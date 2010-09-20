{-# OPTIONS_GHC -fglasgow-exts #-}
-----------------------------------------------------------------------------
-- |
-- Module      :  Numeric.Signal.EEG
-- Copyright   :  (c) Alexander Vivian Hugh McPhail 2010
-- License     :  GPL-style
--
-- Maintainer  :  haskell.vivian.mcphail <at> gmail <dot> com
-- Stability   :  provisional
-- Portability :  
--
-- EEG functions
--
-----------------------------------------------------------------------------

module Numeric.Signal.EEG (
                           loadBDF
                          ) where


import Numeric.Signal.EEG.BDF
