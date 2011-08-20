-- QuasiNewton.hs

{-# LANGUAGE FlexibleContexts #-}
{-# OPTIONS_GHC -Wall #-}

module Optimization.QuasiNewton(dfp, bfgs, dfpHessian, bfgsHessian, quasiNewton, interiorPoint) where

import Optimization.LineSearch
import Numeric.LinearAlgebra
import Casadi.Utils

import System.IO.Unsafe
import Debug.Trace

type HessianApprox a = Matrix a -> Vector a -> Vector a -> Vector a -> Vector a -> Matrix a


dfp :: (Product a, Container Vector a, Ord a, Num (Vector a), Floating a) =>
       Vector a -> (Vector a -> a) -> (Vector a -> Vector a) -> [(Vector a, Matrix a)]
dfp = quasiNewton dfpHessian goldenSectionSearch


bfgs :: (Product a, Container Vector a, Ord a, Num (Vector a), Floating a) =>
       Vector a -> (Vector a -> a) -> (Vector a -> Vector a) -> [(Vector a, Matrix a)]
bfgs = quasiNewton bfgsHessian goldenSectionSearch
--bfgs = quasiNewton bfgsHessian strongWolfe


quasiNewton :: (Product a, Container Vector a, Ord a, Num (Vector a), Floating a) =>
               HessianApprox a -> LineSearch a -> Vector a -> (Vector a -> a) -> (Vector a -> Vector a) -> [(Vector a, Matrix a)]
quasiNewton hessianApprox linesearch x0 f g = iterate (\(x,v) -> 
                                                        if abs (g x) < eps then 
                                                          (x,v) else
                                                          oneQuasiNewton hessianApprox linesearch x v f g) (x0, v0)
  where
    v0 = ident (dim x0)
    eps = 0.0001


oneQuasiNewton :: (Product a, Container Vector a, Ord a, Num (Vector a), Floating a) =>
                  HessianApprox a -> LineSearch a-> Vector a -> Matrix a -> (Vector a -> a) -> (Vector a -> Vector a) -> (Vector a, Matrix a)
oneQuasiNewton hessianApprox linesearch xk vk f g = (xkp1, vkp1)
  where
    vkp1 = hessianApprox vk xk xkp1 gk gkp1
    gk   = g xk
    gkp1 = g xkp1

    pk = -( vk <> gk )
    f' alpha = f $ xk + (scale alpha pk)
    g' alpha = (g $ xk + (scale alpha pk)) <.> pk
    h' _ = pk <.> (vk <> pk) -- constant hessian approximation
    (alphakp1, _) = head $ drop 200 $ linesearch f' g' h' (0,10)
    xkp1 = xk + (scale alphakp1 pk)


bfgsHessian :: (Product a, Container Vector a, Num (Vector a)) => HessianApprox a
bfgsHessian vk xk xkp1 gk gkp1 = vkp1
  where
    vkp1 = (eye - (scale den (outer sk yk))) <> vk <> (eye - (scale den (outer yk sk)))
           + (scale den (outer sk sk))

    eye = ident $ cols vk
    den = 1/(sk <.> yk)

    yk = gkp1 - gk
    sk = xkp1 - xk
    

dfpHessian :: (Product a, Container Vector a, Num (Vector a)) => HessianApprox a
dfpHessian vk xk xkp1 gk gkp1 = vkp1
  where
    vkp1 = vk - ak + bk

    ak = scale (1/(yk <.> (vk <> yk))) (outer (vk <> yk) (yk <> vk))
    bk = scale (1/(sk <.> yk)) (outer sk sk)

    yk = gkp1 - gk
    sk = xkp1 - xk


-- assume f/g and constraints all operate on n-dim vector
-- interiorPoint :: (Floating b, Product a, Container Vector a, Ord a, Num (Vector a), Floating a) =>
--                   Integer -> Vector a -> (Vector a -> a) -> (Vector a -> Vector a)
--               -> [[b] -> b] -> [(Vector a, Matrix a)]
interiorPoint n x0 f g constraints = iterate (\(xs, mu) -> trace ("mu: "++show mu) $ oneInteriorPoint (last xs) (0.9 * mu)) ([x0], 1.0)
  where

    fplus  :: (Num a) => (b -> a) -> (b -> a) -> (b -> a)
    ftimes :: (Num a) => (b -> a) -> (b -> a) -> (b -> a)
    fplus  f g b = (f b) + (g b)
    ftimes f g b = (f b) * (g b)

    barrier  = \v -> sum $ map (negate . log . ($ v)) $ constraints
    (barrierf, barrierg, _) =  unsafePerformIO $ getDerivs barrier n


    oneInteriorPoint x0 mu = 
      let (res, _) = unzip $ bfgs x0 
                     (f `fplus` (const mu `ftimes` barrierf))
                     (g `fplus` ((const (buildVector (fromIntegral n) (\_ -> mu)) `ftimes` barrierg)))
      in (take 22 res, mu)



--    oneInteriorPoint x0 mu = let (res, _) = unzip $ bfgs x0 (f) (g) in (res !! 222, mu)



