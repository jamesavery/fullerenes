module Main where

import Seeds
import Expansion
import qualified Data.IntMap.Strict as IM
import Data.List (sort, nub)

---------------------------------------------------------------------
-- Test infrastructure
---------------------------------------------------------------------

data TestResult = Pass String | Fail String String
    deriving (Show)

passed :: TestResult -> Bool
passed (Pass _) = True
passed _        = False

runTests :: [TestResult] -> IO ()
runTests results = do
    let total = length results
        passes = length (filter passed results)
        fails  = filter (not . passed) results
    mapM_ (\r -> case r of
        Pass name -> putStrLn $ "  PASS: " ++ name
        Fail name msg -> putStrLn $ "  FAIL: " ++ name ++ " -- " ++ msg
        ) results
    putStrLn $ "\n" ++ show passes ++ "/" ++ show total ++ " tests passed."
    if null fails
        then putStrLn "All tests passed!"
        else do
            putStrLn $ show (length fails) ++ " FAILURES:"
            mapM_ (\(Fail n m) -> putStrLn $ "  " ++ n ++ ": " ++ m) fails

---------------------------------------------------------------------
-- Validation
---------------------------------------------------------------------

-- | Check that a DualGraph passes all structural invariants
checkValid :: String -> DualGraph -> TestResult
checkValid name g = case validateDualGraph g of
    Right () -> Pass (name ++ " validates")
    Left err -> Fail (name ++ " validates") err

-- | Check round-trip: expand then reduce should give back the original graph
checkRoundTrip :: String -> DualGraph -> Expansion -> TestResult
checkRoundTrip name g exp =
    let (g', pi) = applyExpansion exp g
        g''      = applyReduction exp pi g'
    in if g'' == g
       then Pass (name ++ " round-trip")
       else Fail (name ++ " round-trip")
            ("Graphs differ after round-trip.\n  Original nv=" ++ show (numVertices g)
            ++ "\n  Expanded nv=" ++ show (numVertices g')
            ++ "\n  Reduced nv=" ++ show (numVertices g''))

-- | Check that expansion produces a valid graph with correct vertex count
checkExpansionValid :: String -> DualGraph -> Expansion -> TestResult
checkExpansionValid name g exp =
    let (g', _) = applyExpansion exp g
    in case validateDualGraph g' of
        Right () -> Pass (name ++ " expansion valid")
        Left err -> Fail (name ++ " expansion valid") err

-- | Check expanded graph has correct vertex count
checkExpansionSize :: String -> DualGraph -> Expansion -> Int -> TestResult
checkExpansionSize name g exp expectedNv =
    let (g', _) = applyExpansion exp g
        actualNv = numVertices g'
    in if actualNv == expectedNv
       then Pass (name ++ " size=" ++ show expectedNv)
       else Fail (name ++ " size=" ++ show expectedNv)
            ("Expected " ++ show expectedNv ++ " vertices, got " ++ show actualNv)

---------------------------------------------------------------------
-- Tests
---------------------------------------------------------------------

-- | Test graph primitives on C20
testPrimitives :: [TestResult]
testPrimitives =
    [ -- C20 has 12 vertices, all degree-5
      check "c20 has 12 vertices" (numVertices c20 == 12)
    , check "c20 all degree-5" (all (\v -> deg c20 v == 5) [0..11])
    , check "c20 validates" (validateDualGraph c20 == Right ())
    -- nextCW/prevCW are inverses
    , check "nextCW . prevCW = id" $
        all (\u -> all (\v -> prevCW c20 u (nextCW c20 u v) == v) (nbrs c20 u)) [0..11]
    , check "prevCW . nextCW = id" $
        all (\u -> all (\v -> nextCW c20 u (prevCW c20 u v) == v) (nbrs c20 u)) [0..11]
    -- advanceCW 0 = identity
    , check "advanceCW 0 = id" $
        all (\u -> all (\v -> advanceCW c20 u v 0 == v) (nbrs c20 u)) [0..11]
    -- advanceCW deg = identity (full cycle)
    , check "advanceCW deg = id" $
        all (\u -> all (\v -> advanceCW c20 u v (deg c20 u) == v) (nbrs c20 u)) [0..11]
    ]
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test L0 expansions on C20
testL0OnC20 :: [TestResult]
testL0OnC20 =
    let exps = expansionsL0 c20
    in [ check ("c20 has " ++ show (length exps) ++ " L0 expansions")
              (length exps > 0)
       -- C20: all 12 vertices are degree-5. Each vertex has 5 neighbors,
       -- all of which are degree-5. Number of edges between degree-5 pairs:
       -- 30 edges (each vertex has 5 neighbors, 12*5/2 = 30). Each with 2 directions = 60.
       , check "c20 has 60 L0 expansions" (length exps == 60)
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_L0_" ++ show i
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp 14  -- 12 + 2 = 14
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 exps))  -- test first 10 for speed
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test that all L0 round-trips work on C20
testAllL0RoundTrips :: [TestResult]
testAllL0RoundTrips =
    let exps = expansionsL0 c20
    in map (\(i, exp) -> checkRoundTrip ("c20_L0_all_" ++ show i) c20 exp)
           (zip [0..] exps)

-- | Test expansions on C28
testC28 :: [TestResult]
testC28 =
    let l0s = expansionsL0 c28
    in [ check "c28 validates" (validateDualGraph c28 == Right ())
       , check "c28 has 16 vertices" (numVertices c28 == 16)
       , check ("c28 has " ++ show (length l0s) ++ " L0 expansions") True
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c28_L0_" ++ show i
        in [ checkExpansionValid name c28 exp
           , checkExpansionSize name c28 exp 18
           , checkRoundTrip name c28 exp
           ]) (zip [0..] (take 5 l0s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test straight expansions (L_i, i >= 1) on expanded C20
testStraightOnExpandedC20 :: [TestResult]
testStraightOnExpandedC20 =
    -- First expand C20 via L0 to get a 14-vertex graph, then look for L1 expansions
    let exp0 = head (expansionsL0 c20)
        (g14, _) = applyExpansion exp0 c20
    in [ checkValid "c20+L0" g14
       , check "c20+L0 has 14 vertices" (numVertices g14 == 14)
       ]
    ++ let l1s = expansionsL 4 g14
       in [ check ("c20+L0 has " ++ show (length l1s) ++ " L_i expansions") True ]
       ++ concatMap (\(i, exp) ->
            let name = "c20+L0_L" ++ show (expI exp) ++ "_" ++ show i
            in [ checkExpansionValid name g14 exp
               , checkRoundTrip name g14 exp
               ]) (zip [0..] (take 5 l1s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"
    expI (Exp (L i) _ _) = i
    expI _ = -1

-- | Test B_{0,0} on C20
testBentZeroOnC20 :: [TestResult]
testBentZeroOnC20 =
    let b00s = [ Exp (B 0 0) (u, v) d
               | u <- degree5 c20
               , v <- nbrs c20 u
               , d <- [DLeft, DRight]
               , canBentPath c20 (u, v) d 0 0
               ]
    in [ check ("c20 has " ++ show (length b00s) ++ " B00 candidates") True ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_B00_" ++ show i
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp 15  -- 12 + 3 = 15
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 b00s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test general bent expansions on C20
testBentOnC20 :: [TestResult]
testBentOnC20 =
    let bents = expansionsB 4 c20
        -- Filter to just B_{i,j} with i+j > 0
        nonZeroBents = filter (\(Exp (B i j) _ _) -> i + j > 0) bents
    in [ check ("c20 has " ++ show (length nonZeroBents) ++ " B_{i,j} (i+j>0) expansions") True ]
    ++ concatMap (\(i, exp) ->
        let name = "c20_B_" ++ show i
            Exp (B bi bj) _ _ = exp
        in [ checkExpansionValid name c20 exp
           , checkExpansionSize name c20 exp (12 + bi + bj + 3)
           , checkRoundTrip name c20 exp
           ]) (zip [0..] (take 10 nonZeroBents))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test on C30
testC30 :: [TestResult]
testC30 =
    let l0s = expansionsL0 c30
    in [ check "c30 validates" (validateDualGraph c30 == Right ())
       , check "c30 has 17 vertices" (numVertices c30 == 17)
       , check ("c30 has " ++ show (length l0s) ++ " L0 expansions") True
       ]
    ++ concatMap (\(i, exp) ->
        let name = "c30_L0_" ++ show i
        in [ checkExpansionValid name c30 exp
           , checkRoundTrip name c30 exp
           ]) (zip [0..] (take 5 l0s))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test bent expansions on a larger graph (expand C28 to get ~20 vertices)
testBentOnExpanded :: [TestResult]
testBentOnExpanded =
    -- Expand C28 via a couple of L0s to get a larger graph
    let exp0 = head (expansionsL0 c28)
        (g18, _) = applyExpansion exp0 c28
        exp1 = head (expansionsL0 g18)
        (g20, _) = applyExpansion exp1 g18
    in [ checkValid "c28+2L0" g20
       , check ("c28+2L0 has " ++ show (numVertices g20) ++ " vertices") (numVertices g20 == 20)
       ]
    ++ let bents = expansionsB 4 g20
       in [ check ("c28+2L0 has " ++ show (length bents) ++ " bent expansions") True ]
       ++ concatMap (\(i, exp) ->
            let name = "c28+2L0_B_" ++ show i
                Exp (B bi bj) _ _ = exp
            in [ checkExpansionValid name g20 exp
               , checkExpansionSize name g20 exp (numVertices g20 + bi + bj + 3)
               , checkRoundTrip name g20 exp
               ]) (zip [0..] (take 10 bents))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

-- | Test all expansion types comprehensively on an expanded graph
testAllExpansionTypes :: [TestResult]
testAllExpansionTypes =
    let -- Build a larger graph by expanding C20 multiple times
        exp0 = head (expansionsL0 c20)
        (g14, _) = applyExpansion exp0 c20
        exp1 = head (expansionsL0 g14)
        (g16, _) = applyExpansion exp1 g14
        exp2 = head (expansionsL0 g16)
        (g18, _) = applyExpansion exp2 g16
        -- Now test ALL expansion types on g18
        allExps = expansions 4 g18
        l0s   = [e | e@(Exp (L 0) _ _) <- allExps]
        lis   = [e | e@(Exp (L i) _ _) <- allExps, i > 0]
        bents = [e | e@(Exp (B _ _) _ _) <- allExps]
    in [ checkValid "c20+3L0" g18
       , check ("c20+3L0 has " ++ show (numVertices g18) ++ " vertices") (numVertices g18 == 18)
       , check ("All expansions: " ++ show (length allExps)
                ++ " (L0=" ++ show (length l0s)
                ++ " L_i=" ++ show (length lis)
                ++ " B=" ++ show (length bents) ++ ")") True
       ]
    ++ concatMap (\(i, exp) ->
        let Exp k _ _ = exp
            name = "g18_" ++ show k ++ "_" ++ show i
        in [ checkExpansionValid name g18 exp
           , checkRoundTrip name g18 exp
           ]) (zip [0..] (take 20 allExps))
  where
    check name True  = Pass name
    check name False = Fail name "condition is False"

---------------------------------------------------------------------
-- Main
---------------------------------------------------------------------

main :: IO ()
main = do
    putStrLn "=== Graph Primitives ==="
    runTests testPrimitives

    putStrLn "\n=== L0 Expansions on C20 ==="
    runTests testL0OnC20

    putStrLn "\n=== All L0 Round-Trips on C20 ==="
    runTests testAllL0RoundTrips

    putStrLn "\n=== C28 Tests ==="
    runTests testC28

    putStrLn "\n=== Straight Expansions on Expanded C20 ==="
    runTests testStraightOnExpandedC20

    putStrLn "\n=== B_{0,0} on C20 ==="
    runTests testBentZeroOnC20

    putStrLn "\n=== General Bent on C20 ==="
    runTests testBentOnC20

    putStrLn "\n=== C30 Tests ==="
    runTests testC30

    putStrLn "\n=== Bent Expansions on Expanded C28 ==="
    runTests testBentOnExpanded

    putStrLn "\n=== All Expansion Types on Expanded C20 ==="
    runTests testAllExpansionTypes
