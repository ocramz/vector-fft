import Test.Tasty
import Test.Tasty.QuickCheck hiding ((.&.))
import Test.QuickCheck.Instances.Vector ()

import Data.Bits
import Data.Complex
import qualified Data.Vector.Unboxed as V
import Data.Vector.FFT

expi :: Double -> Complex Double
expi y = cos y :+ sin y

dft :: V.Vector (Complex Double) -> V.Vector (Complex Double)
dft v = V.generate len (\k -> V.sum $ V.imap (\n x -> x * expi (-2 * pi * fromIntegral k * fromIntegral n / fromIntegral len)) v)
  where
    len = V.length v

idft :: V.Vector (Complex Double) -> V.Vector (Complex Double)
idft v = V.generate len (\k -> (V.sum $ V.imap (\n x -> x * expi (2 * pi * fromIntegral k * fromIntegral n / fromIntegral len)) v) / fromIntegral len)
  where
    len = V.length v

(~=~) :: V.Vector (Complex Double) -> V.Vector (Complex Double) -> Property
v1 ~=~ v2 = counterexample (show v1 ++ (if res then " == " else " /= ") ++ show v2) res
  where
    epsilon = 1e-10
    res = V.all (\(x, y) -> magnitude (x - y) < epsilon) (V.zip v1 v2)

-- | Next power of 2
nextPow2 :: Int -> Int
nextPow2 n
  -- `n .&. (n - 1)` clears the lowest set bit
  | n .&. (n - 1) == 0 = n
  | otherwise = 1 `shiftL` (log2 n + 1)

log2 :: Int -> Int
log2 n = finiteBitSize (0 :: Int) - 1 - countLeadingZeros n

pad :: V.Vector (Complex Double) -> V.Vector (Complex Double)
pad v = V.generate (nextPow2 len) (\i -> if i < len then v V.! i else 0)
  where
    len = V.length v

main :: IO ()
main = defaultMain $ localOption (QuickCheckTests 1000) $ testGroup "vector-fft"
    [ testProperty "ifft . fft = id" $ \v -> ifft (fft v) ~=~ v
    , testProperty "fft . ifft = id" $ \v -> fft (ifft v) ~=~ v
    , testProperty "fft = dft" $ \v -> fft v ~=~ dft (pad v)
    , testProperty "ifft = idft" $ \v -> ifft v ~=~ idft (pad v)
    ]
