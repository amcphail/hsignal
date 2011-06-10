{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
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
                Filterable(..),
                freqz,
                pwelch,
                hilbert
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
type PF = Ptr Float

-----------------------------------------------------------------------------

class Convolvable a where
    -- | convolve two containers, output is the size of the second argument, no zero padding
    convolve :: a -> a -> a

-----------------------------------------------------------------------------

class (Storable a, Container Vector a, Num (Vector a)
      , Convert a, Floating (Vector a), RealElement a
      , Num a)
       => Filterable a where
    -- | conver from Vector Double
    fromDouble :: Vector Double -> Vector a
--       b ~ ComplexOf a, Container Vector b, Convert b) => Filterable a where
    -- | filter a signal
    filter_ :: Vector a -- ^ zero coefficients
            -> Vector a -- ^ pole coefficients
            -> Vector a -- ^ input signal
            -> Vector a -- ^ output signal
    -- | coefficients of a Hamming window
    hamming_ :: Int           -- ^ length
           -> Vector a -- ^ the Hamming coeffficents
    -- | the complex power : real $ v * (conj v)
    complex_power_ :: Vector (Complex Double) -- ^ input
                  -> Vector a                       -- ^ output
     -- | resample, take one sample every n samples in the original
    downsample_ :: Int -> Vector a -> Vector a
    -- | the difference between consecutive elements of a vector
    deriv_ :: Vector a -> Vector a
    -- | unwrap the phase of signal (input expected to be within (-pi,pi)
    unwrap_ :: Vector a -> Vector a
    -- | evaluate a real coefficient polynomial for complex arguments
    polyEval_ :: Vector a           -- ^ the real coefficients
         -> Vector (Complex Double) -- ^ the points at which to be evaluated
         -> Vector (Complex Double) -- ^ the values
    -- | the cross covariance of two signals
    cross_covariance_ :: Int            -- ^ maximum delay
                     -> Vector a       -- ^ time series
                     -> Vector a       -- ^ time series
                     -> (a,a,Vector a) -- ^ (sd_x,sd_y,cross_covariance)
    -- | the cumulative sum of a signal
    cumulative_sum_ :: Vector a    -- ^ time series
                    -> Vector a    -- ^ result

-----------------------------------------------------------------------------

instance Convolvable (Vector Double) where
    convolve x y = fst $ fromComplex $ F.ifft $ (F.fft (complex x) * F.fft (complex y))
--    convolve = convolve_vector_double

convolve_vector_double c a = unsafePerformIO $ do
                             r <- createVector (dim a)
                             app3 signal_vector_double_convolve vec c vec a vec r "signalDoubleConvolve"
                             return r

foreign import ccall "signal-aux.h vector_double_convolve" signal_vector_double_convolve :: CInt -> PD -> CInt -> PD -> CInt -> PD -> IO CInt

instance Convolvable (Vector Float) where
    convolve x y = single $ fst $ fromComplex $ F.ifft $ (F.fft (complex $ double x) * F.fft (complex $ double y))
--    convolve = convolve_vector_double

convolve_vector_float c a = unsafePerformIO $ do
                             r <- createVector (dim a)
                             app3 signal_vector_float_convolve vec c vec a vec r "signalFloatConvolve"
                             return r

foreign import ccall "signal-aux.h vector_float_convolve" signal_vector_float_convolve :: CInt -> PF -> CInt -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

instance Convolvable (Vector (Complex Double)) where
    convolve x y = F.ifft $ (F.fft x * F.fft y)
--    convolve = convolve_vector_complex

convolve_vector_complex c a = unsafePerformIO $ do
                              r <- createVector (dim a)
                              app3 signal_vector_complex_convolve vec c vec a vec r "signalComplexConvolve"
                              return r

foreign import ccall "signal-aux.h vector_complex_convolve" signal_vector_complex_convolve :: CInt -> PC -> CInt -> PC -> CInt -> PC -> IO CInt

instance Convolvable (Vector (Complex Float)) where
    convolve x y = single $ F.ifft $ (F.fft (double x) * F.fft (double y))

-----------------------------------------------------------------------------

instance Filterable Double where
    fromDouble = id
    filter_ = filterD
    hamming_ = hammingD
    complex_power_ = complex_powerD
    downsample_ = downsampleD
    deriv_ = derivD
    unwrap_ = unwrapD
    polyEval_ = polyEval
    cross_covariance_ = crossCovarianceD
    cumulative_sum_ = cumSumD

instance Filterable Float where
    fromDouble = single
    filter_ = filterF
    hamming_ = hammingF
    complex_power_ = complex_powerF
    downsample_ = downsampleF
    deriv_ = derivF
    unwrap_ = unwrapF
    polyEval_ c = polyEval (double c)
    cross_covariance_ = crossCovarianceF
    cumulative_sum_ = cumSumF

-----------------------------------------------------------------------------

-- | filters the signal
filterD :: Vector Double -- ^ zero coefficients
       -> Vector Double -- ^ pole coefficients
       -> Vector Double -- ^ input signal
       -> Vector Double -- ^ output signal
filterD l k v = unsafePerformIO $ do
               r <- createVector (dim v)
               app4 signal_filter_double vec l vec k vec v vec r "signalFilter"
               return r

foreign import ccall "signal-aux.h filter_double" signal_filter_double :: CInt -> PD -> CInt -> PD -> CInt -> PD -> CInt -> PD -> IO CInt

-- | filters the signal
filterF :: Vector Float -- ^ zero coefficients
       -> Vector Float -- ^ pole coefficients
       -> Vector Float -- ^ input signal
       -> Vector Float -- ^ output signal
filterF l k v = unsafePerformIO $ do
               r <- createVector (dim v)
               app4 signal_filter_float vec l vec k vec v vec r "signalFilter"
               return r

foreign import ccall "signal-aux.h filter_float" signal_filter_float :: CInt -> PF -> CInt -> PF -> CInt -> PF -> CInt -> PF -> IO CInt

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
hammingD :: Int           -- ^ length
        -> Vector Double -- ^ the Hamming coeffficents
hammingD l 
    | l == 1          = constant 1.0 1
    | otherwise       = unsafePerformIO $ do
                        r <- createVector l
                        app1 signal_hamming_double vec r "Hamming"
                        return r

foreign import ccall "signal-aux.h hamming_double" signal_hamming_double :: CInt -> PD -> IO CInt

-- | coefficients of a Hamming window
hammingF :: Int           -- ^ length
        -> Vector Float -- ^ the Hamming coeffficents
hammingF l 
    | l == 1          = constant 1.0 1
    | otherwise       = unsafePerformIO $ do
                        r <- createVector l
                        app1 signal_hamming_float vec r "Hamming"
                        return r

foreign import ccall "signal-aux.h hamming_float" signal_hamming_float :: CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

-- | determine the frequency response of a filter
{-freqz :: (Filterable a, Storable a, Container Vector a, Convert a, RealElement a,
         DoubleOf a ~ DoubleOf (RealOf b), RealElement c, c ~ DoubleOf a, c ~ DoubleOf (RealOf b),
         b ~ Complex a, b ~ ComplexOf a, Convert b, Container Vector b,
         Container Vector c, Convert c, b ~ ComplexOf a, b ~ ComplexOf c)
        ⇒ Vector a     -- ^ zero coefficients
      -> Vector a       -- ^ pole coefficients
      -> Vector a       -- ^ points (between 0 and 2*pi)
      -> Vector a       -- ^ response
-}
freqz :: (Filterable a, Complex Double ~ ComplexOf (DoubleOf a)
        ,Filterable (DoubleOf a)) => 
        Vector a       -- ^ zero coefficients
      -> Vector a       -- ^ pole coefficients
      -> Vector a       -- ^ points (between 0 and 2*pi)
      -> Vector a       -- ^ response
freqz b a w = let k = max (dim b) (dim a)
                  hb = polyEval_ (postpad b k) (exp (scale (0 :+ 1) ((complex $ double w))))
                  ha = polyEval_ (postpad a k) (exp (scale (0 :+ 1) ((complex $ double w))))
              in complex_power_ (hb / ha)

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
complex_powerD :: Vector (Complex Double) -- ^ input
              -> Vector Double           -- ^ output
complex_powerD v = unsafePerformIO $ do
                  r <- createVector (dim v)
                  app2 signal_complex_power_double vec v vec r "complex_power"
                  return r

foreign import ccall "signal-aux.h complex_power_double" signal_complex_power_double :: CInt -> PC -> CInt -> PD -> IO CInt

-- | the complex power : real $ v * (conj v)
complex_powerF :: Vector (Complex Double) -- ^ input
              -> Vector Float             -- ^ output
complex_powerF v = unsafePerformIO $ do
                  r <- createVector (dim v)
                  app2 signal_complex_power_float vec v vec r "complex_power"
                  return r

foreign import ccall "signal-aux.h complex_power_float" signal_complex_power_float :: CInt -> PC -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

-- | resample, take one sample every n samples in the original
downsampleD :: Int -> Vector Double -> Vector Double
downsampleD n v = unsafePerformIO $ do
               r <- createVector (dim v `div` n)
               app2 (signal_downsample_double $ fromIntegral n) vec v vec r "downsample"
               return r

foreign import ccall "signal-aux.h downsample_double" signal_downsample_double :: CInt -> CInt -> PD -> CInt -> PD -> IO CInt

-- | resample, take one sample every n samples in the original
downsampleF :: Int -> Vector Float -> Vector Float
downsampleF n v = unsafePerformIO $ do
               r <- createVector (dim v `div` n)
               app2 (signal_downsample_float $ fromIntegral n) vec v vec r "downsample"
               return r

foreign import ccall "signal-aux.h downsample_float" signal_downsample_float :: CInt -> CInt -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

-- | the difference between consecutive elements of a vector
derivD :: Vector Double -> Vector Double
derivD v = unsafePerformIO $ do
          r <- createVector (dim v - 1)
          app2 (signal_diff_double) vec v vec r "diff"
          return r

foreign import ccall "signal-aux.h vector_diff_double" signal_diff_double :: CInt -> PD -> CInt -> PD -> IO CInt

-- | the difference between consecutive elements of a vector
derivF :: Vector Float -> Vector Float
derivF v = unsafePerformIO $ do
          r <- createVector (dim v - 1)
          app2 (signal_diff_float) vec v vec r "diff"
          return r

foreign import ccall "signal-aux.h vector_diff_float" signal_diff_float :: CInt -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

-- | unwrap the phase of signal (input expected to be within (-pi,pi)
unwrapD :: Vector Double -> Vector Double
unwrapD v = unsafePerformIO $ do
           r <- createVector $ dim v
           app2 signal_unwrap_double vec v vec r "unwrap"
           return r

foreign import ccall "signal-aux.h unwrap_double" signal_unwrap_double :: CInt -> PD -> CInt -> PD -> IO CInt

-- | unwrap the phase of signal (input expected to be within (-pi,pi)
unwrapF :: Vector Float -> Vector Float
unwrapF v = unsafePerformIO $ do
           r <- createVector $ dim v
           app2 signal_unwrap_float vec v vec r "unwrap"
           return r

foreign import ccall "signal-aux.h unwrap_float" signal_unwrap_float :: CInt -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

-- | compute the cross covariance of two signals
crossCovarianceD :: Int -> Vector Double -> Vector Double -> (Double,Double,Vector Double)
crossCovarianceD l x y = unsafePerformIO $ do
                           r <- createVector (2*l)
                           alloca $ \sx -> 
                               alloca $ \sy -> do
                                 app3 (signal_cross_covariance_double (fromIntegral l) sx sy) vec x vec y vec r "cross_covariance"
                                 sx' <- peek sx
                                 sy' <- peek sy
                                 return (sx',sy',r)

foreign import ccall "signal-aux.h cross_covariance_double" 
        signal_cross_covariance_double :: CInt -> PD -> PD -> CInt -> PD -> CInt
                                       -> PD -> CInt -> PD -> IO CInt

-- | compute the cross covariance of two signals
crossCovarianceF :: Int -> Vector Float -> Vector Float -> (Float,Float,Vector Float)
crossCovarianceF l x y = unsafePerformIO $ do
                           r <- createVector (2*l)
                           alloca $ \sx -> 
                               alloca $ \sy -> do
                                 app3 (signal_cross_covariance_float (fromIntegral l) sx sy) vec x vec y vec r "cross_covariance"
                                 sx' <- peek sx
                                 sy' <- peek sy
                                 return (sx',sy',r)

foreign import ccall "signal-aux.h cross_covariance_float" 
        signal_cross_covariance_float :: CInt -> PF -> PF -> CInt -> PF -> CInt
                                       -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------

cumSumD :: Vector Double -> Vector Double
cumSumD v = unsafePerformIO $ do
              r <- createVector (dim v)
              app2 signal_cum_sum_double vec v vec r "cumSumD"
              return r

cumSumF :: Vector Float -> Vector Float
cumSumF v = unsafePerformIO $ do
              r <- createVector (dim v)
              app2 signal_cum_sum_float vec v vec r "cumSumF"
              return r

foreign import ccall "signal-aux.h cum_sum_double"
        signal_cum_sum_double :: CInt -> PD -> CInt -> PD -> IO CInt
foreign import ccall "signal-aux.h cum_sum_float"
        signal_cum_sum_float :: CInt -> PF -> CInt -> PF -> IO CInt

-----------------------------------------------------------------------------



