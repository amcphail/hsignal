{-# OPTIONS_GHC -fglasgow-exts #-}
{-# OPTIONS_GHC -XUndecidableInstances #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Signal.Multichannel
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  uses Concurrency
--
-- Signal processing functions, multichannel datatype
--
-- link with '-threaded' and run with +RTS -Nn, where n is the number of CPUs
--
-----------------------------------------------------------------------------

module Numeric.Signal.Multichannel (
                       Multichannel,readMultichannel,writeMultichannel,
                       createMultichannel,
                       sampling_rate,precision,channels,samples,
                       detrended,filtered,
                       getChannel,getChannels,
                       toMatrix,
                       mapConcurrently,
                       detrend,filter,
                       slice,
                       mi_phase
                ) where

-----------------------------------------------------------------------------

import qualified Numeric.Signal as S

--import Complex

import qualified Data.Array.IArray as I
import Data.Ix

import Control.Concurrent
--import Control.Concurrent.MVar

import System.IO.Unsafe(unsafePerformIO)

--import qualified Data.List as L

import Data.Binary

import Foreign.Storable

import Numeric.LinearAlgebra

--import qualified Numeric.GSL.Fourier as F

import qualified Numeric.GSL.Histogram as H
import qualified Numeric.GSL.Histogram2D as H2

import qualified Numeric.Statistics.Information as SI

import Prelude hiding(filter)

import Control.Monad(replicateM)


{-
-------------------------------------------------------------------

instance (Binary a, Storable a) => Binary (Vector a) where
    put v = do
            let d = dim v
            put d
            mapM_ (\i -> put $ v @> i) [0..(d-1)]
    get = do
          d <- get
          xs <- replicateM d get
          return $ fromList xs

-------------------------------------------------------------------
-}

-----------------------------------------------------------------------------

-- | data type with multiple channels
data Multichannel a = MC {
                          _sampling_rate :: Int             -- ^ sampling rate
                          , _precision   :: Int             -- ^ bits of precision
                          , _channels    :: Int             -- ^ number of channels
                          , _length      :: Int             -- ^ length in samples
                          , _detrended   :: Bool            -- ^ was the data detrended?
                          , _filtered    :: Maybe (Int,Int) -- ^ if filtered the passband
                          , _data        :: I.Array Int (Vector a) -- ^ data
                         }

-----------------------------------------------------------------------------

instance (Binary a, Storable a, 
          Ord a, RealFrac a,
          Container Vector a,
         Product a) => Binary (Multichannel a) where
    put (MC s p c l de f d) = do
                              put s
                              put p
                              put c
                              put l
                              put de
                              put f
                              put $! fmap convert d
        where convert v = let (mi,ma) = (minElement v,maxElement v)
                              v' = mapVector (\x -> round $ (x - mi)/(ma - mi) * (fromIntegral (maxBound :: Word32))) v
                          in (mi,ma,(v' :: Vector Word32)) 



    get = do
          s <- get
          p <- get
          c <- get
          l <- get
          de <- get
          f <- get
          d <- (get :: Get (I.Array Int (a,a,Vector Word32)))
          return $! (MC s p c l de f (seq d (fmap convert) d))
              where convert (mi,ma,v) = mapVector (\x -> ((fromIntegral x) :: a) / (fromIntegral (maxBound :: Word32)) * (ma - mi) + mi) v

-----------------------------------------------------------------------------

readMultichannel :: (Binary a, Storable a, 
                     Ord a, RealFrac a,
                     Container Vector a,
                    Product a) => FilePath -> IO (Multichannel a)
readMultichannel = decodeFile

writeMultichannel :: (Binary a, Storable a, 
                      Ord a, RealFrac a,
                      Container Vector a,
                     Product a) => FilePath -> Multichannel a -> IO ()
writeMultichannel = encodeFile

-----------------------------------------------------------------------------

-- | create a multichannel data type
createMultichannel :: Storable a 
                   => Int               -- ^ sampling rate
                   -> Int               -- ^ bits of precision
                   -> [Vector a]        -- ^ data
                   -> Multichannel a    -- ^ datatype
createMultichannel s p d = let c = length d
                 in MC s p c (dim $ head d) False Nothing (I.listArray (1,c) d)

-- | the sampling rate
sampling_rate :: Multichannel a -> Int
sampling_rate = _sampling_rate

-- | the bits of precision
precision :: Multichannel a -> Int
precision = _precision

-- | the number of channels
channels :: Multichannel a -> Int
channels = _channels

-- | the length, in samples
samples :: Multichannel a -> Int
samples = _length

-- | extract one channel
getChannel :: Int -> Multichannel a -> Vector a
getChannel c d = (_data d) I.! c

-- | extract all channels
getChannels :: Multichannel a -> I.Array Int (Vector a)
getChannels d = _data d

-- | convert the data to a matrix with channels as rows
toMatrix :: Element a => Multichannel a -> Matrix a
toMatrix = fromRows . I.elems . _data

-- | was the data detrended?
detrended :: Multichannel a -> Bool
detrended = _detrended

-- | was the data filtered?
filtered :: Multichannel a -> Maybe (Int,Int)
filtered = _filtered

-----------------------------------------------------------------------------

-- | map a function executed concurrently
mapArrayConcurrently :: Ix i => (a -> b)    -- ^ function to map
                     -> I.Array i a         -- ^ input
                     -> I.Array i b         -- ^ output
mapArrayConcurrently f d = unsafePerformIO $ do
                                             let b = I.bounds d
                                             results <- replicateM (rangeSize b) newEmptyMVar
                                             mapM_ (forkIO . applyFunction f) $ zip results (I.assocs d)
                                             vectors <- mapM takeMVar results
                                             return $ I.array (I.bounds d) vectors
    where applyFunction f' (m,(j,e)) = putMVar m (j,f' e)

{-
-- | map a function executed concurrently
mapListConcurrently :: (a -> b)             -- ^ function to map
                    -> [a]                  -- ^ input
                    -> [b]                  -- ^ output
mapListConcurrently f d = unsafePerformIO $ do
                                            results <- replicateM (length d) newEmptyMVar
                                            mapM_ (forkIO . applyFunction f) zip results d
                                            mapM takeMVar results
    where applyFunction f' (m,e) = putMVar m (f' e)
-}

-- | map a function executed concurrently
mapConcurrently :: Storable b 
                => (Vector a -> Vector b)      -- ^ the function to be mapped 
                -> Multichannel a              -- ^ input data
                -> Multichannel b              -- ^ output data
mapConcurrently f (MC sr p c _ de fi d) = let d' = mapArrayConcurrently f d
                                          in MC sr p c (dim $ d' I.! 1) de fi d'

-- | map a function
mapMC :: Storable b 
      => (Vector a -> Vector b)                -- ^ the function to be mapped 
      -> Multichannel a                        -- ^ input data
      -> Multichannel b                        -- ^ output data
mapMC f (MC sr p c _ de fi d) = let d' = fmap f d
                                in MC sr p c (dim $ d' I.! 1) de fi d'
                                    
-----------------------------------------------------------------------------

-- | detrend the data with a specified window size
detrend :: Int -> Multichannel Double -> Multichannel Double
detrend w m = let m' = mapConcurrently (S.detrend w) m
              in m' { _detrended = True }


-- | filter the data with the given passband
filter :: (S.Filterable a, Double ~ DoubleOf a, Container Vector (Complex a), Convert (Complex a)) ⇒ 
         (Int,Int) -> Multichannel a -> Multichannel a
filter pb m = let m' = mapConcurrently (S.broadband_filter (_sampling_rate m) pb) m
              in m' { _filtered = Just pb }

-----------------------------------------------------------------------------

-- | extract a slice of the data
slice :: Storable a 
      => Int                 -- ^ starting sample number
      -> Int                 -- ^ length
      -> Multichannel a 
      -> Multichannel a
slice j w m = let m' = mapConcurrently (subVector j w) m
              in m' { _length = w }

-----------------------------------------------------------------------------

-- | calculate the mutual information of the phase between pairs of channels (fills upper half of matrix)
mi_phase :: (S.Filterable a, Double ~ DoubleOf a) ⇒
           Multichannel a      -- ^ input data
         -> Matrix Double
mi_phase m = let d = fmap double $ _data m
                 histarray = mapArrayConcurrently (H.fromLimits 128 (-pi,pi)) d
                 c = channels m
                 pairs = I.array ((1,1),(c,c)) $ map (\(a,b) -> ((a,b),((a,b),d I.! a,d I.! b))) (range ((1,1),(c,c)))
                 hist2array = mapArrayConcurrently (\(j,x,y) -> (j,H2.addVector (H2.emptyLimits 128 128 (-pi,pi) (-pi,pi)) x y)) pairs
                 mi = mapArrayConcurrently (doMI histarray d) hist2array
             in fromArray2D mi
    where doMI histarray d ((x,y),h2) 
              | x < y     = SI.mutual_information h2 (histarray I.! x) (histarray I.! y) (d I.! x,d I.! y)
              | otherwise = 0

-----------------------------------------------------------------------------
