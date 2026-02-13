module Main where

import Seeds
import qualified Data.IntMap.Strict as IM

main :: IO ()
main = do
    putStrLn "Testing seed graph decoding and validation..."
    putStrLn ""

    -- General seeds
    putStrLn "=== General seeds ==="
    testSeed "C20" c20 12
    testSeed "C28" c28 16
    testSeed "C30" c30 17

    -- IPR seeds with valid boundary
    putStrLn "\n=== IPR seeds (valid boundary) ==="
    let iprExpectedNv = [46, 44, 38, 42, 42, 47, 40, 32, 42, 42]
    mapM_ (\(i, (g, nv)) -> testSeed ("IPR[" ++ show i ++ "]") g nv)
          (zip [0..] (zip iprSeeds iprExpectedNv))

    -- IPR seeds with invalid boundary
    putStrLn "\n=== IPR seeds (invalid boundary) ==="
    let iprInvExpectedNv = [51, 43, 56, 44, 46, 50, 45, 58, 44, 43, 43, 47, 48, 48,
                            44, 44, 44, 44, 44, 44, 44, 44, 43, 43, 43, 43, 43,
                            42, 42, 41, 40, 39, 37, 42, 41, 46]
    mapM_ (\(i, (g, nv)) -> testSeed ("IPR_inv[" ++ show i ++ "]") g nv)
          (zip [0..] (zip iprSeedsInvalidBoundary iprInvExpectedNv))

    putStrLn "\n====================================="
    putStrLn "All 49 seeds decoded and validated successfully!"

testSeed :: String -> DualGraph -> Int -> IO ()
testSeed name g expectedNv = do
    -- Check vertex count
    let nv = numVertices g
    check (name ++ " vertex count") (nv == expectedNv)
        ("expected " ++ show expectedNv ++ ", got " ++ show nv)

    -- Check edge count (triangulation: E = 3V - 6)
    let adj = neighbours g
        totalDeg = sum [length (adj IM.! v) | v <- [0..nv-1]]
        e = totalDeg `div` 2
        expectedE = 3 * nv - 6
    check (name ++ " edge count") (e == expectedE)
        ("expected " ++ show expectedE ++ ", got " ++ show e)

    -- Check degree-5 count
    let n5 = length (degree5 g)
    check (name ++ " degree-5 count") (n5 == 12)
        ("expected 12, got " ++ show n5)

    -- Full validation
    case validateDualGraph g of
        Left err -> fail $ name ++ " validation failed: " ++ err
        Right () -> putStrLn $ "  " ++ name ++ ": OK (V=" ++ show nv ++
                               ", E=" ++ show e ++ ", deg5=" ++ show n5 ++ ")"

check :: String -> Bool -> String -> IO ()
check label True  _   = return ()
check label False msg = fail $ label ++ ": " ++ msg
