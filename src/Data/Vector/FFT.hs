{-# language BangPatterns #-}
{-# options_ghc -Wno-unused-imports #-}
module Data.Vector.FFT (
  fft, ifft
  -- * Useful results
  , crossCorrelation
  ) where

import Control.Monad (when)
import Control.Monad.Primitive (PrimMonad(..))

import Control.Monad.ST (runST)
import Data.Bits (countLeadingZeros, finiteBitSize, shiftR, shiftL, (.&.), (.|.))
import Data.Bool (Bool, otherwise)
import Data.Complex (Complex(..), conjugate)
import Data.Foldable (forM_)

import Data.Vector.Unboxed as V (Vector, Unbox, map, zipWith, length, unsafeFreeze, (!))
import qualified Data.Vector.Unboxed.Mutable as VM (MVector, read, write, new, length)
import qualified Data.Vector.Generic as VG (Vector(..), copy)

import Prelude hiding (read)


-- | (Circular) cross-correlation of two vectors
--
-- Defined via the FFT and IFFT for computational efficiency
--
-- NB the source vectors should have matching length for meaningful results
crossCorrelation :: Vector (Complex Double)
                 -> Vector (Complex Double)
                 -> Vector (Complex Double)
crossCorrelation v1 v2 = ifft $ (cmap conjugate v1hat) `prod` v2hat
  where
    prod = V.zipWith (*)
    v1hat = fft v1
    v2hat = fft v2



-- | Radix-2 decimation-in-time fast Fourier Transform.
--
--   The given array (and therefore the output as well) is zero-padded to the next power of two if necessary.
fft :: Vector (Complex Double) -> Vector (Complex Double)
fft arr = runST $ do
  marr <- copyPadded arr
  mfft marr
  V.unsafeFreeze marr
{-# inlinable [1] fft #-}

-- | Inverse fast Fourier transform.
--
--   The given array (and therefore the output as well) is zero-padded to the next power of two if necessary.
ifft :: Vector (Complex Double) -> Vector (Complex Double)
ifft arr = do
  let lenComplex = intToComplexDouble (nextPow2 (V.length arr))
  cmap ((/ lenComplex) . conjugate) . fft . cmap conjugate $ arr
{-# inlinable [1] ifft #-}




-- | Copy the source vector into a zero-padded mutable one
copyPadded :: (PrimMonad m, Num a, Unbox a) =>
              Vector a -> m (VM.MVector (PrimState m) a)
copyPadded arr = do
  let
    len = V.length arr
    l2 = nextPow2 len
  marr <- VM.new l2
  forM_ [0 .. l2 - 1] $ \i -> do
    let x | i < len = arr V.! i
          | otherwise = 0
    VM.write marr i x
  pure marr
{-# inline copyPadded #-}





-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two,
--   though this property is not checked.
mfft :: (PrimMonad m) => VM.MVector (PrimState m) (Complex Double) -> m ()
mfft mut = do {
    let len = VM.length mut
  ; let bitReverse !i !j = do {
          ; if i == len - 1
              then stage 0 1
              else do {
                  when (i < j) $ swap mut i j
                ; let inner k l = if k <= l
                        then inner (k `shiftR` 1) (l - k)
                        else bitReverse (i + 1) (l + k)
                ; inner (len `shiftR` 1) j
              }
        }
        stage l l1 = if l == (log2 len)
          then pure ()
          else do {
              let !l2 = l1 `shiftL` 1
                  !e = (negate twoPi) / (intToDouble l2)
                  flight j !a = if j == l1
                    then stage (l + 1) l2
                    else do {
                        let butterfly i = if i >= len
                              then flight (j + 1) (a + e)
                              else do {
                                  let i1 = i + l1
                                ; xi1 :+ yi1 <- VM.read mut i1
                                ; let !co = cos a
                                      !si = sin a
                                      d = (co * xi1 - si * yi1) :+ (si * xi1 + co * yi1)
                                ; ci <- VM.read mut i
                                ; VM.write mut i1 (ci - d)
                                ; VM.write mut i (ci + d)
                                ; butterfly (i + l2)
                              }
                      ; butterfly j
                    }
            ; flight 0 0
         }
  ; when (len > 0) $ bitReverse 0 0
}

-- | Next power of 2
nextPow2 :: Int -> Int
nextPow2 n
  -- `n .&. (n - 1)` clears the lowest set bit
  | n .&. (n - 1) == 0 = n
  | otherwise = 1 `shiftL` (log2 n + 1)


log2 :: Int -> Int
log2 n = finiteBitSize (0 :: Int) - 1 - countLeadingZeros n


{-# inline swap #-}
swap :: (PrimMonad m, Unbox a) =>
        VM.MVector (PrimState m) a -> Int -> Int -> m ()
swap mut i j = do
  atI <- VM.read mut i
  atJ <- VM.read mut j
  VM.write mut i atJ
  VM.write mut j atI

twoPi :: Double
{-# inline twoPi #-}
twoPi = 6.283185307179586

intToDouble :: Int -> Double
{-# inline intToDouble #-}
intToDouble = fromIntegral

intToComplexDouble :: Int -> Complex Double
{-# inline intToComplexDouble #-}
intToComplexDouble = fromIntegral


{-# inline cmap #-}
cmap :: (Floating a, Unbox a) => (Complex a -> Complex a) -> V.Vector (Complex a) -> V.Vector (Complex a)
cmap = V.map
