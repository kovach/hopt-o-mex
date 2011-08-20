-- optTest.hs

{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleInstances, UndecidableInstances #-}

module Main where

import Casadi.Utils
import Optimization.QuasiNewton
import Numeric.LinearAlgebra
import Data.Packed.Vector

import Graphics.Gnuplot.Simple

rosenbrock :: Floating a => [a] -> a
rosenbrock [x0,x1] = (1-x0)*(1-x0) + 100*(x1-x0*x0)*(x1-x0*x0)
rosenbrock _ = error "wrong number of inputs to rosenbrock"

quadratic [x,y] = x^2+y^2

class (Show a, Read a) => Lol a
instance (Show a, Read a) => Lol a

shovel :: (Lol a) => a -> String
shovel a = show a

zero [x0,x1] = x0
zero _ = error "wrong"



main :: IO ()
main = do
  (rosenbrockF, rosenbrockG, _) <- getDerivs rosenbrock 2
  (quadf, quadg, _) <- getDerivs quadratic 2
  (zf, zg, _) <- getDerivs zero 2
  let x0 = fromList [-0.9,1.1] :: Vector Double
--      results = (take 20 $ dfp x0 rosenbrockF rosenbrockG)
--      results = (take 20 $ bfgs x0 rosenbrockF rosenbrockG)
--      results = take 40 $ interiorPoint 2 x0 rosenbrockF rosenbrockG [\[x,y] -> (0.5-x)]
      results = take 40 $ interiorPoint 2 x0 quadf quadg [\[x,y] -> y-1]
      path = concatMap f $ tail results
        where
          f (vecs,_) = map (\ [x,y] -> (x,y)) $ map toList vecs

--  mapM_ (\(_,v) -> print $ inv v) results
--  putStrLn "hi"
--  print $ zf $ fromList [3,3]

  mapM_ (\(x,_) -> print x) results
  
  plotList [] path