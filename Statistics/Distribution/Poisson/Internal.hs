-- |
-- Module    : Statistics.Distribution.Poisson.Internal
-- Copyright : (c) 2011 Bryan O'Sullivan
-- License   : BSD3
--
-- Maintainer  : bos@serpentine.com
-- Stability   : experimental
-- Portability : portable
--
-- Internal code for the Poisson distribution.

module Statistics.Distribution.Poisson.Internal
    (
      probability,
      logProbability
    ) where

import Numeric.MathFunctions.Constants (m_sqrt_2_pi, m_tiny, m_neg_inf)
import Numeric.SpecFunctions           (logGammaFN, stirlingError)
import Numeric.SpecFunctions.Extra     (bd0)

-- | An unchecked, non-integer-valued version of Loader's saddle point
-- algorithm.
probability :: Double -> Double -> Double
probability 0      0     = 1
probability 0      1     = 0
probability lambda x
  | isInfinite lambda    = 0
  | x < 0                = 0
  | x <= lambda * m_tiny = exp (-lambda)
  | lambda < x * m_tiny  = exp (-lambda + x * log lambda - logGammaFN (x+1))
  | otherwise            = exp (-(stirlingError x) - bd0 x lambda) /
                           (m_sqrt_2_pi * sqrt x)
{-# INLINE probability #-}

logProbability :: Double -> Double -> Double
logProbability 0      0     = 0
logProbability 0      1     = m_neg_inf
logProbability lambda x
  | isInfinite lambda    = m_neg_inf
  | x < 0                = m_neg_inf
  | x <= lambda * m_tiny = -lambda
  | lambda < x * m_tiny  = -lambda + x * log lambda - logGammaFN (x+1)
  | otherwise            = (-(stirlingError x) - bd0 x lambda) - 0.5 * (log (2*pi*x))
{-# INLINE logProbability #-}
