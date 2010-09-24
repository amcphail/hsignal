{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Signal.Internal
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  uses FFI
--
-- low-level interface
--
-----------------------------------------------------------------------------

module Numeric.Signal.Internal (
                Convolvable(..),
                hamming,
                filter,
                freqz,
                pwelch,
                hilbert,
                complex_power,
                downsample,
                deriv,
                unwrap
                ) where

import Data.Packed.Development(createVector,vec,app1,app2,app3,app4)
--import Data.Packed.Vector
import Numeric.LinearAlgebra

--import Numeric.LinearAlgebra.Algorithms
--import Numeric.LinearAlgebra.Linear

import qualified Numeric.GSL.Fourier as F
import Foreign
--import Data.Complex
import Foreign.C.Types

import Prelude hiding(filter)

-----------------------------------------------------------------------------

type PD = Ptr Double                            
type PC = Ptr (Complex Double)                  

-----------------------------------------------------------------------------

class Convolvable a where
    -- | convolve two containers, output is the size of the second argument, no zero padding
    convolve :: a -> a -> a

-----------------------------------------------------------------------------

instance Convolvable (Vector Double) where
    convolve x y = fst $ fromComplex $ F.ifft $ (F.fft (complex x) * F.fft (complex y))
--    convolve = convolve_vector_double

convolve_vector_double c a = unsafePerformIO $ do
                             r <- createVector (dim a)
                             app3 signal_vector_double_convolve vec c vec a vec r "signalDoubleConvolve"
                             return r

foreign import ccall "signal-aux.h vector_double_convolve" signal_vector_double_convolve :: CInt -> PD -> CInt -> PD -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

instance Convolvable (Vector (Complex Double)) where
    convolve x y = F.ifft $ (F.fft x * F.fft y)
--    convolve = convolve_vector_complex

convolve_vector_complex c a = unsafePerformIO $ do
                              r <- createVector (dim a)
                              app3 signal_vector_complex_convolve vec c vec a vec r "signalComplexConvolve"
                              return r

foreign import ccall "signal-aux.h vector_complex_convolve" signal_vector_complex_convolve :: CInt -> PC -> CInt -> PC -> CInt -> PC -> IO CInt

-----------------------------------------------------------------------------

-- | filters the signal
filter :: Vector Double -- ^ zero coefficients
       -> Vector Double -- ^ pole coefficients
       -> Vector Double -- ^ input signal
       -> Vector Double -- ^ output signal
filter l k v = unsafePerformIO $ do
               r <- createVector (dim v)
               app4 signal_filter vec l vec k vec v vec r "signalFilter"
               return r

foreign import ccall "signal-aux.h filter" signal_filter :: CInt -> PD -> CInt -> PD -> CInt -> PD -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | Hilbert transform with original vector as real value, transformed as imaginary
hilbert :: Vector Double -> Vector (Complex Double)
hilbert v = unsafePerformIO $ do
            let r = complex v
            -- could use (complex v) to make a complex vector in haskell rather than C
            app1 signal_hilbert vec r "hilbert"
            return r

foreign import ccall "signal-aux.h hilbert" signal_hilbert :: CInt -> PC -> IO CInt

-----------------------------------------------------------------------------

-- | Welch (1967) power spectrum density using periodogram/FFT method
pwelch :: Int            -- ^ window size (multiple of 2)
       -> Vector Double  -- ^ input signal
       -> Vector Double  -- ^ power density  
pwelch w v = unsafePerformIO $ do
             let r = constant 0.0 ((w `div` 2) + 1)
             app2 (signal_pwelch $ fromIntegral w) vec (complex v) vec r "pwelch"
             return r

foreign import ccall "signal-aux.h pwelch" signal_pwelch :: CInt -> CInt -> PC -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | coefficients of a Hamming window
hamming :: Int           -- ^ length
        -> Vector Double -- ^ the Hamming coeffficents
hamming l 
    | l == 1          = constant 1.0 1
    | otherwise       = unsafePerformIO $ do
                        r <- createVector l
                        app1 signal_hamming vec r "Hamming"
                        return r

foreign import ccall "signal-aux.h hamming" signal_hamming :: CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | determine the frequency response of a filter
freqz :: Vector Double     -- ^ zero coefficients
      -> Vector Double     -- ^ pole coefficients
      -> Vector Double     -- ^ points (between 0 and 2*pi)
      -> Vector Double     -- ^ response
freqz b a w = let k = max (dim b) (dim a)
                  hb = polyEval (postpad b k) (exp (scale (0 :+ 1) ((complex w) :: Vector (Complex Double))))
                  ha = polyEval (postpad a k) (exp (scale (0 :+ 1) ((complex w) :: Vector (Complex Double))))
              in complex_power (hb / ha)

postpad v n = let d = dim v
              in if d < n then join [v,(constant 0.0 (n-d))]
              else v

-----------------------------------------------------------------------------

-- | evaluate a real coefficient polynomial for complex arguments
polyEval :: Vector Double           -- ^ the real coefficients
         -> Vector (Complex Double) -- ^ the points at which to be evaluated
         -> Vector (Complex Double) -- ^ the values
polyEval c z = unsafePerformIO $ do
               r <- createVector (dim z)
               app3 signal_real_poly_complex_eval vec c vec z vec r "polyEval"
               return r

foreign import ccall "signal-aux.h real_poly_complex_eval" signal_real_poly_complex_eval :: CInt -> PD -> CInt -> PC -> CInt -> PC -> IO CInt

-----------------------------------------------------------------------------

-- | the complex power : real $ v * (conj v)
complex_power :: Vector (Complex Double) -- ^ input
              -> Vector Double           -- ^ output
complex_power v = unsafePerformIO $ do
                  r <- createVector (dim v)
                  app2 signal_complex_power vec v vec r "complex_power"
                  return r

foreign import ccall "signal-aux.h complex_power" signal_complex_power :: CInt -> PC -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | resample, take one sample every n samples in the original
downsample :: Int -> Vector Double -> Vector Double
downsample n v = unsafePerformIO $ do
               r <- createVector (dim v `div` n)
               app2 (signal_downsample $ fromIntegral n) vec v vec r "downsample"
               return r

foreign import ccall "signal-aux.h downsample" signal_downsample :: CInt -> CInt -> PD -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | the difference between consecutive elements of a vector
deriv :: Vector Double -> Vector Double
deriv v = unsafePerformIO $ do
          r <- createVector (dim v - 1)
          app2 (signal_deriv) vec v vec r "diff"
          return r

foreign import ccall "signal-aux.h vector_deriv" signal_deriv :: CInt -> PD -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------

-- | unwrap the phase of signal (input expected to be within (-pi,pi)
unwrap :: Vector Double -> Vector Double
unwrap v = unsafePerformIO $ do
           r <- createVector $ dim v
           app2 signal_unwrap vec v vec r "unwrap"
           return r

foreign import ccall "signal-aux.h unwrap" signal_unwrap :: CInt -> PD -> CInt -> PD -> IO CInt

-----------------------------------------------------------------------------



