{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Signal
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  uses FFI
--
-- Signal processing functions
--
-----------------------------------------------------------------------------

module Numeric.Signal (
                       -- * Filtering
                       hamming,
                       pwelch,
                       fir,standard_fir,broadband_fir,
                       freqzF,freqzN,
                       filter,broadband_filter,
                       -- * Analytic Signal
                       analytic_signal,analytic_power,analytic_phase,
                       unwrap,
                       -- * Preprocessing
                       detrend,
                       downsample,resize,
                       -- * Utility functions
                       deriv
                ) where

-----------------------------------------------------------------------------

import qualified Numeric.Signal.Internal as S

import Numeric.GSL.Fitting.Linear
 
import Data.Complex

import qualified Data.List as L

--import Data.Packed.Vector
--import Data.Packed(Container(..))
import Numeric.LinearAlgebra

import qualified Numeric.GSL.Fourier as F

import Prelude hiding(filter)

-----------------------------------------------------------------------------

-- | filters the signal
filter :: Vector Double -- ^ zero coefficients
       -> Vector Double -- ^ pole coefficients
       -> Int           -- ^ sampling rate
       -> Vector Double -- ^ input signal
       -> Vector Double -- ^ output signal
filter b a s v = let len = dim v
                     w = min s len
                     start = (negate . fromList . reverse . toList . subVector 0 w) v
                     finish = (negate . fromList . reverse . toList . subVector (len-w) w) v
                     v' = join [start,v,finish]
                 in subVector s len $ S.filter b a v'

-----------------------------------------------------------------------------
                     
-- | coefficients of a Hamming window
hamming :: Int           -- ^ length
        -> Vector Double -- ^ the Hamming coeffficents
hamming = S.hamming

-----------------------------------------------------------------------------

-- | Welch (1967) power spectrum density using periodogram/FFT method
pwelch :: Int            -- ^ sampling rate
       -> Int            -- ^ window size
       -> Vector Double  -- ^ input signal
       -> (Vector Double,Vector Double)  -- ^ (frequency index,power density)  
pwelch s w v = let w' = max s w -- make window at least sampling rate
                   r  = S.pwelch w' v
                   sd = (fromIntegral s)/2
                   -- scale for sampling rate
                   r' = scale (recip sd) r
                   f  = linspace ((w `div` 2) + 1) (0,sd)
               in (f,r')

-----------------------------------------------------------------------------

-- | a broadband FIR
broadband_fir :: Int           -- ^ sampling rate
              -> (Int,Int)     -- ^ (lower,upper) frequency cutoff
              -> Vector Double -- ^ filter coefficients   
broadband_fir s (l,h) = let o = 501
                            ny = (fromIntegral s) / 2.0
                            fl = (fromIntegral l) / ny
                            fh = (fromIntegral h) / ny
                            f = [0, fl*0.95, fl, fh, fh*1.05, 1]
                            m = [0,0,1,1,0,0]
                            be = zip f m
                        in standard_fir o be

-- | a broadband filter
broadband_filter :: Int           -- ^ sampling rate
                 -> (Int,Int)     -- ^ (lower,upper) frequency cutoff
                 -> Vector Double -- ^ input signal
                 -> Vector Double -- ^ output signal
broadband_filter s f v = let b = broadband_fir s f
                         in filter b (scalar 1.0) s v
                                
-----------------------------------------------------------------------------

-- | standard FIR filter
-- |   FIR filter with grid a power of 2 greater than the order, ramp = grid/16, hamming window
standard_fir :: Int -> [(Double,Double)] -> Vector Double
standard_fir o be = let grid  = calc_grid o
                        trans = grid `div` 16
                    in fir o be grid trans $ S.hamming (o+1)

calc_grid :: Int -> Int
calc_grid o = let next_power = ceiling (((log $ fromIntegral o) :: Double) / (log 2.0)) :: Int
              in floor $ 2.0 ** ((fromIntegral next_power) :: Double)


-- | produce an FIR filter
fir :: Int               -- ^ order (one less than the length of the filter)
    -> [(Double,Double)] -- ^ band edge frequency, nondecreasing, [0, f1, ..., f(n-1), 1]
                         -- ^ band edge magnitude
    -> Int               -- ^ grid spacing
    -> Int               -- ^ transition width
    -> Vector Double     -- ^ smoothing window (size is order + 1)
    -> Vector Double     -- ^ the filter coefficients
fir o be gn tn w = let mid = o `div` 2
                       (f,m) = unzip be
                       f' = diff (((fromIntegral gn) :: Double)/((fromIntegral tn) :: Double)/2.0) f
                       m' = interpolate f m f'
                       grid = interpolate f' m' $ map (\x -> (fromIntegral x)/(fromIntegral gn)) [0..(gn-1)]
                       grid' = map (\x -> x :+ 0) grid
                       (b,_) = fromComplex $ F.ifft $ fromList $ grid' ++ (reverse (drop 1 grid'))
                       b' = join [subVector ((dim b)-mid-1) (mid+1) b, subVector 1 (mid+1) b] 
                   in b' * w

floor_zero x
    | x < 0.0   = 0.0
    | otherwise = x

ceil_one x
    | x > 1.0   = 1.0
    | otherwise = x

diff :: Double -> [Double] -> [Double]
diff _ []  = []
diff _ [x] = [x]
diff inc (x1:x2:xs)
     | x1 == x2     = (floor_zero $ x1-inc):x1:(ceil_one $ x1+inc):(diff inc (L.filter (/= x2) xs))
     | otherwise    = x1:(diff inc (x2:xs))

interpolate :: [Double] -> [Double] -> [Double] -> [Double]
interpolate _ _ []      = []
interpolate x y (xp:xs) = if xp == 1.0 
                             then ((interpolate'' ((length x)-1) x y xp):(interpolate x y xs))
                             else ((interpolate' x y xp):(interpolate x y xs))

interpolate' :: [Double] -> [Double] -> Double -> Double
interpolate' x y xp = let Just j = L.findIndex (> xp) x
                      in (interpolate'' j x y xp)

interpolate'' :: Int -> [Double] -> [Double] -> Double -> Double
interpolate'' j x y xp = let x0 = x !! (j-1)
                             y0 = y !! (j-1)
                             x1 = x !! j
                             y1 = y !! j
                         in y0 + (xp - x0) * ((y1 - y0)/(x1-x0))

-----------------------------------------------------------------------------

-- | determine the frequency response of a filter, given a vector of frequencies
freqzF :: Vector Double     -- ^ zero coefficients
       -> Vector Double     -- ^ pole coefficients
       -> Int               -- ^ sampling rate   
       -> Vector Double     -- ^ frequencies
       -> Vector Double     -- ^ frequency response
freqzF b a s f = S.freqz b a ((2*pi/(fromIntegral s)) * f)

-- | determine the frequency response of a filter, given a number of points and sampling rate
freqzN :: Vector Double     -- ^ zero coefficients
       -> Vector Double     -- ^ pole coefficients
       -> Int               -- ^ sampling rate
       -> Int               -- ^ number of points
       -> (Vector Double,Vector Double)     -- ^ (frequencies,response)
freqzN b a s n = let w' = linspace n (0,((fromIntegral n)-1)/(fromIntegral (2*n)))
                     r = S.freqz b a ((2*pi)*w')
                     in ((fromIntegral s)*w',r)
                     
-----------------------------------------------------------------------------

-- | an analytic signal is the original signal with Hilbert-transformed signal as imaginary component
analytic_signal :: Vector Double -> Vector (Complex Double)
analytic_signal = S.hilbert

-- | the power (amplitude^2 = v * (conj c)) of an analytic signal
analytic_power :: Vector (Complex Double) -> Vector Double
analytic_power = S.complex_power

-- | the phase of an analytic signal
analytic_phase :: Vector (Complex Double) -> Vector Double
analytic_phase = uncurry arctan2 . fromComplex

-----------------------------------------------------------------------------

-- | remove a linear trend from data
detrend :: Int             -- ^ window size
        -> Vector Double   -- ^ data to be detrended
        -> Vector Double   -- ^ detrended data
detrend w v = let windows = dim v `div` w
                  ws = takesV ((replicate windows w) ++ [dim v - (windows * w)]) v
                  ds = map detrend' ws
                  windows' = (dim v - (w `div` 2)) `div` w
                  ws' = takesV (((w `div` 2):(replicate windows' w)) ++ [dim v - (w `div` 2) - (windows' * w)]) v
                  ds' = map detrend' ws'
              in (join ds + join ds') / 2 
    where detrend' x = let ln = dim x
                           t = linspace ln (1.0,fromIntegral ln)
                           (c0,c1,_,_,_,_) = linear t x
                       in x - (scale c1 t + scalar c0)

-----------------------------------------------------------------------------

-- | take one sample from every n samples in the original
downsample :: Int -> Vector Double -> Vector Double
downsample = S.downsample


-- | resize the vector to length n by resampling
resize :: Int -> Vector Double -> Vector Double
resize n v = downsample (dim v `div` n) v

-----------------------------------------------------------------------------

-- | the difference between consecutive elements of a vector
deriv :: Vector Double -> Vector Double
deriv = S.deriv

-----------------------------------------------------------------------------

-- | unwrap the phase of signal (input expected to be within (-pi,pi)
unwrap :: Vector Double -> Vector Double
unwrap = S.unwrap

-----------------------------------------------------------------------------
