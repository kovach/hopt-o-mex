-- DdpCasadi.hs

{-# OPTIONS_GHC -Wall #-}
{-# LANGUAGE FlexibleContexts, RankNTypes #-}

module DdpCasadi
       (
         prepareQFunction,
         evalQFunction,
         prepareDdp,
         BacksweepOutput
       ) where

import Casadi
import Hom(State, Action, Cost, Ode, Quad(..), evalQuad)
import Numeric.LinearAlgebra
import Data.List(mapAccumL)
import System.IO.Unsafe(unsafePerformIO)

type BacksweepOutput a = (Quad a, [[a]], [a])

-- prepare casadi SXFunction
prepareQFunction :: Integral a => a -> a -> Cost SX -> Ode SX -> IO SXFunction
prepareQFunction nx nu costFunction dode = do
  let nx' = fromIntegral nx
      nu' = fromIntegral nu
  x   <- sxMatrixCreateSymbolic "x"   (nx',1)
  u   <- sxMatrixCreateSymbolic "u"   (nu',1)
  vxx <- sxMatrixCreateSymbolic "vxx" (nx',nx')
  vx  <- sxMatrixCreateSymbolic "vx"  (nx',1)
  v0  <- sxMatrixCreateSymbolic "v0"  (1,1)
  x0  <- sxMatrixCreateSymbolic "x0"  (nx',1)
      
  let -- reshape SXMatrices into lists
      x'   = sxMatrixToList x
      u'   = sxMatrixToList u
      vxx' = sxMatrixToLists vxx
      vx'  = sxMatrixToList vx
      v0'  = head (sxMatrixToList v0)
      x0'  = sxMatrixToList x0
      quad = Quad vxx' vx' v0' x0'
      
      -- the qFunction
      q = sxMatrixFromList [(costFunction x' u') + (evalQuad quad (dode x' u'))]
      qFun = sxFunctionCreate [x,u,vxx,vx,v0,x0] [q]
      
      -- the quadratic expansion
      qx  = sxFunctionGradientAt qFun 0
      qu  = sxFunctionGradientAt qFun 1
      qxx = sxFunctionHessianAt qFun (0,0)
      qxu = sxFunctionHessianAt qFun (0,1)
      quu = sxFunctionHessianAt qFun (1,1)

  return $ sxFunctionCreate [x,u,vxx,vx,v0,x0] [q, qx, qu, qxx, qxu, quu]

-- evaluate
--evalQFunction :: (Real a) => SXFunction -> (State a, Action a, Quad a) -> IO (a,[a],[a],[[a]],[[a]],[[a]])
evalQFunction :: SXFunction -> (State Double, Action Double, Quad Double) -> IO ([[Double]],[[Double]],[[Double]],[[Double]],[[Double]],[[Double]])
evalQFunction qExpansion (x, u, Quad vxx vx v0 x0) = do
  [q,qx,qu,qxx,qxu,quu] <- sxFunctionEvaluate qExpansion [[x],[u],vxx,[vx],[[v0]],[x0]]
  return (q,qx,qu,qxx,qxu,quu)


----------------------- convenience function -----------------------
-- prepare ddp  
prepareDdp :: Integral a => (forall b. Floating b => Cost b) -> (forall b. Floating b => Ode b) -> a -> a
       -> IO ([State Double] -> [Action Double] -> [([State Double], [Action Double], [BacksweepOutput Double])])
prepareDdp cost dode nx nu = do
  qFun <- prepareQFunction nx nu cost dode
  return $ ddp qFun dode
  
  
-- iterate ddp'
ddp :: SXFunction -> Ode Double -> [State Double] -> [Action Double] -> [([State Double], [Action Double], [BacksweepOutput Double])]
ddp cost dode xTraj0 uTraj0 = iterate f x0
  where
    f (xTraj, uTraj, _) = ddp' cost dode xTraj uTraj
    x0 = ddp' cost dode xTraj0 uTraj0


-- just one backsweep then forwardsweep
ddp' :: SXFunction -> Ode Double -> [State Double] -> [Action Double] -> ([State Double], [Action Double], [BacksweepOutput Double])
ddp' qFunction dode xTraj0 uTraj0 = (xTraj, uTraj, backsweepTrajectory)
  where
    backsweepTrajectory :: [BacksweepOutput Double]
    backsweepTrajectory = backSweep qFunction xTraj0 uTraj0

    forwardsweepTrajectory :: [(State Double, Action Double)]
    forwardsweepTrajectory = forwardSweep dode (head xTraj0) backsweepTrajectory
    (xTraj, uTraj) = unzip forwardsweepTrajectory


-------------------------- forward sweep -----------------------------
forwardSweep :: (Floating a) => Ode a -> State a -> [BacksweepOutput a]
                -> [(State a, Action a)]
forwardSweep dode x0 backsweepTrajectory = snd $ mapAccumL (fSweep dode) x0 backsweepTrajectory


fSweep :: Floating a => Ode a -> State a -> BacksweepOutput a -> (State a, (State a, Action a))
fSweep dode x (Quad _ _ _ x0, feedbackGain, uOpenLoop) = (xNext, (x, u))
  where
    u = mimoController x x0 feedbackGain uOpenLoop
    xNext = dode x u

mimoController :: Floating a => State a -> State a -> [[a]] -> Action a -> Action a
mimoController x x0 feedbackMatrix uOpenLoop = u
  where
    u = zipWith (+) uOpenLoop uClosedLoop
    uClosedLoop = mvMult feedbackMatrix deltaX
      where
        deltaX = zipWith (-) x x0
        mvMult m v = map (\a -> dotLists a v) m
          where
            dotLists xa xb = sum (zipWith (*) xa xb)


-------------------------- backward sweep --------------------------
backSweep :: SXFunction -> [State Double] -> [Action Double] -> [BacksweepOutput Double]
backSweep qFunction xTraj0 uTraj0 = foldr (\x acc -> backSweep' qFunction x acc) [] (zip xTraj0 uTraj0)

backSweep' :: SXFunction -> (State Double, Action Double) -> [(Quad Double, [[Double]], [Double])] -> [BacksweepOutput Double]
-- end step - Vnext is 0 so q fcn is cost function
backSweep' qFunction (x,u) [] = [backPropagate qFunction x u (Quad vxx vx v0 x0)]
  where
    x0 = vx
    vxx = replicate (length x) vx
    vx = replicate (length x) v0
    v0 = 0
-- all non-end steps
backSweep' qFunction (x,u) acc@((v,_,_):_) = (backPropagate qFunction x u v):acc


-- convert lists to matrix/vectors, apply the q function back propogation math, convert back to lists
backPropagate :: SXFunction -> State Double -> Action Double -> Quad Double -> BacksweepOutput Double
backPropagate qFunction x u nextValue = (Quad vxx vx v0 x, feedbackGains, openLoopControl)
  where
    -- q functions to lists to matrices/vectors
    (q0', qx', qu', qxx', qxu', quu') = unsafePerformIO $ do (evalQFunction qFunction (x,u,nextValue))
    qxx = fromLists qxx'
    quu = fromLists quu'
    qxu = fromLists qxu'
    qx = fromList $ concat qx'
    qu = fromList $ concat qu'
    q0 = head (head q0')
    
    -- value function
    vxx = toLists $ qxx - (qxu <> (inv quu) <> (trans qxu))
    vx = toList $ qx - ( qu <> (inv quu) <> (trans qxu))
    v0 = q0 - (qu `dot` ((inv quu) <> qu))
    
    -- feedback gain
    feedbackGains = toLists $ - ( inv quu) <> (trans qxu);
    
    -- open loop control
    openLoopControl = zipWith (+) u $ toList $ - (inv quu) <> qu;
