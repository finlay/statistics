{-# LANGUAGE CPP, FlexibleContexts, Rank2Types #-}
-- |
-- Module    : Statistics.Function
-- Copyright : (c) 2009, 2010, 2011 Bryan O'Sullivan
-- License   : BSD3
--
-- Maintainer  : bos@serpentine.com
-- Stability   : experimental
-- Portability : portable
--
-- Useful functions.

module Statistics.Function
    (
    -- * Scanning
      minMax
    -- * Sorting
    , sort
    , sortBy
    , partialSort
    -- * Indexing
    , indexed
    , indices
    -- * Bit twiddling
    , nextHighestPowerOfTwo
    -- * Comparison
    , within
    ) where

#include "MachDeps.h"

import Data.Bits ((.|.), shiftR)
import qualified Data.Vector.Algorithms.Intro as I
import qualified Data.Vector.Generic as G
import Statistics.Function.Comparison (within)

-- | Sort a vector.
sort :: (Ord e, G.Vector v e) => v e -> v e
sort = G.modify I.sort
{-# INLINE sort #-}

-- | Sort a vector using a custom ordering.
sortBy :: (G.Vector v e) => I.Comparison e -> v e -> v e
sortBy f = G.modify $ I.sortBy f
{-# INLINE sortBy #-}

-- | Partially sort a vector, such that the least /k/ elements will be
-- at the front.
partialSort :: (G.Vector v e, Ord e) =>
               Int -- ^ The number /k/ of least elements.
            -> v e
            -> v e
partialSort k = G.modify (`I.partialSort` k)
{-# INLINE partialSort #-}

-- | Return the indices of a vector.
indices :: (G.Vector v a, G.Vector v Int) => v a -> v Int
indices a = G.enumFromTo 0 (G.length a - 1)
{-# INLINE indices #-}

-- | Zip a vector with its indices.
indexed :: (G.Vector v e, G.Vector v Int, G.Vector v (Int,e)) => v e -> v (Int,e)
indexed a = G.zip (indices a) a
{-# INLINE indexed #-}

data MM = MM {-# UNPACK #-} !Double {-# UNPACK #-} !Double

-- | Compute the minimum and maximum of a vector in one pass.
minMax :: (G.Vector v Double) => v Double -> (Double, Double)
minMax = fini . G.foldl' go (MM (1/0) (-1/0))
  where
    go (MM lo hi) k = MM (min lo k) (max hi k)
    fini (MM lo hi) = (lo, hi)
{-# INLINE minMax #-}

-- | Efficiently compute the next highest power of two for a
-- non-negative integer.  If the given value is already a power of
-- two, it is returned unchanged.  If negative, zero is returned.
nextHighestPowerOfTwo :: Int -> Int
nextHighestPowerOfTwo n = o + 1
    where m = n - 1
          o = m
              .|. (m `shiftR` 1)
              .|. (m `shiftR` 2)
              .|. (m `shiftR` 4)
              .|. (m `shiftR` 8)
              .|. (m `shiftR` 16)
#if WORD_SIZE_IN_BITS == 64              
              .|. (m `shiftR` 32)
#endif                
