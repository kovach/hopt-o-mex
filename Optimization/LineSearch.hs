-- LineSearch.hs

{-# OPTIONS_GHC -Wall #-}

module Optimization.LineSearch
       (
         goldenSectionSearch
       , strongWolfe
       , LineSearch
       ) where

import Numeric.GSL.Differentiation
import Debug.Trace


d f = fst . (derivCentral 0.001 f)

tau :: Floating a => a
tau = 2/(1+sqrt(5))

type LineSearch a = (a -> a) -> (a -> a) -> (a -> a) -> (a, a) -> [(a,a)]

goldenSectionSearch :: (Ord a, Floating a) => LineSearch a
goldenSectionSearch f _ _ (x0,x3) = gss f (x0,x1,x2,x3)
  where
    x1 = x0 + (x3-x0)*(1-tau)
    x2 = x0 + (x3-x0)*tau

gss :: (Ord a, Floating a) => (a -> a) -> (a, a, a, a) -> [(a,a)]
gss f (x0, x1, x2, x3) = xs
  where
    x1' = x0 + (x2-x0)*(1-tau)
    x2' = x1 + (x3-x1)*tau
    xs
      | (f x1) < (f x2) = (x1, f x1):(gss f (x0,x1',x1,x2))
      | otherwise       = (x2, f x2):(gss f (x1,x2,x2',x3))

--strongWolfe :: (Ord a, Floating a) => LineSearch a
strongWolfe phi gradient _ (alpha1, alphamax) = trace "LOLOL" $ repeat (loop 1 alpha1 alpha0, 0)
  where
    alpha0 = 0
--    alpha1 = 1.0
--    alphamax = 22
    mu1 = 0.0001 
    mu2 = 0.9 
    d _ = gradient
    dphi0 = d phi 0
    loop i alpha oldalpha = 
      if phi alpha > phi 0 + mu1 * alpha * dphi0 ||
         (phi alpha > phi oldalpha && i > 1)
      then 
        zoom oldalpha alpha
      else
        let dphi = d phi alpha in
        {-trace (show (phi alpha) ++ "/" ++ show (dphi0) ++ "\n"
                   ++ show dphi ++ "/" ++ show (-mu2*dphi0)) $-}
        if abs dphi <= -mu2 * dphi0 then alpha else
          if dphi >= 0 then zoom oldalpha alpha else
{-            trace (show alpha) $-} loop (i+1) ((alpha+alphamax) / 2) alpha

    zoom a b = 
      let alpha = (a+b)/2
          dphi = d phi alpha
      in
        if phi alpha > phi 0 + mu1 * alpha * dphi0 ||
           phi alpha > phi a then zoom a alpha
        else
          if abs dphi <= -mu2 * dphi0 then alpha else
            if dphi * (b - a) >= 0 then zoom a alpha else zoom alpha b

