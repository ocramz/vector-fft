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
import qualified Data.Vector.Unboxed.Mutable as VM (MVector, generate, read, write, swap, new, length)
import qualified Data.Vector.Generic as VG (Vector(..), copy)

import Prelude hiding (read)


-- | (Circular) cross-correlation of two vectors.
--
-- Defined via the FFT and IFFT for computational efficiency.
--
-- NB: the source vectors should have matching length for meaningful results.
crossCorrelation :: (RealFloat a, Unbox a) => Vector (Complex a) -> Vector (Complex a) -> Vector (Complex a)
crossCorrelation v1 v2 = ifft $ (cmap conjugate v1hat) `prod` v2hat
  where
    prod = V.zipWith (*)
    v1hat = fft v1
    v2hat = fft v2
{-# specialize crossCorrelation :: Vector (Complex Double) -> Vector (Complex Double) -> Vector (Complex Double) #-}
{-# specialize crossCorrelation :: Vector (Complex Float) -> Vector (Complex Float) -> Vector (Complex Float) #-}

-- | Radix-2 decimation-in-time fast Fourier Transform.
--
--   The given array (and therefore the output as well) is zero-padded to the next power of two if necessary.
fft :: (RealFloat a, Unbox a) => Vector (Complex a) -> Vector (Complex a)
fft arr = runST $ do
  marr <- copyPadded arr
  mfft marr
  V.unsafeFreeze marr
{-# specialize fft :: Vector (Complex Double) -> Vector (Complex Double) #-}
{-# specialize fft :: Vector (Complex Float) -> Vector (Complex Float) #-}

-- | Inverse fast Fourier transform.
--
--   The given array (and therefore the output as well) is zero-padded to the next power of two if necessary.
ifft :: (RealFloat a, Unbox a) => Vector (Complex a) -> Vector (Complex a)
ifft arr = do
  let lenComplex = intToComplex (nextPow2 (V.length arr))
  cmap ((/ lenComplex) . conjugate) . fft . cmap conjugate $ arr
{-# specialize ifft :: Vector (Complex Double) -> Vector (Complex Double) #-}
{-# specialize ifft :: Vector (Complex Float) -> Vector (Complex Float) #-}

-- | Copy the source vector into a zero-padded mutable one
copyPadded :: (PrimMonad m, Num a, Unbox a) => Vector a -> m (VM.MVector (PrimState m) a)
copyPadded arr = do
  let len = V.length arr
  VM.generate (nextPow2 len) $ \i -> if i < len then arr V.! i else 0
{-# inline copyPadded #-}

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two,
--   though this property is not checked.
mfft :: (PrimMonad m, RealFloat a, Unbox a) => VM.MVector (PrimState m) (Complex a) -> m ()
mfft mut = do
    let
      len = VM.length mut

      -- bit reversal permutation
      -- i is the original index
      -- j is the bit-reversed index
      reverseIndices !i !j
        | i == len - 1 = stage 0 1
        | otherwise = do
            when (i < j) $ VM.swap mut i j
            let inner k l
                  | k <= l = inner (k `shiftR` 1) (l - k)
                  | otherwise = reverseIndices (i + 1) (l + k)
            inner (len `shiftR` 1) j

      stage l l1
        | l == (log2 len) = pure ()
        | otherwise = do
            let !l2 = l1 `shiftL` 1
                !e = (-2 * pi) / (intToRealFloat l2)
                flight j !a
                  | j == l1 = stage (l + 1) l2
                  | otherwise = do
                      let butterfly i
                            | i >= len = flight (j + 1) (a + e)
                            | otherwise = do
                                let i1 = i + l1
                                xi1 :+ yi1 <- VM.read mut i1
                                let !co = cos a
                                    !si = sin a
                                    d = (co * xi1 - si * yi1) :+ (si * xi1 + co * yi1)
                                ci <- VM.read mut i
                                VM.write mut i1 (ci - d)
                                VM.write mut i (ci + d)
                                butterfly (i + l2)
                      butterfly j
            flight 0 0

    when (len > 0) $ reverseIndices 0 0

-- | Next power of 2
nextPow2 :: Int -> Int
nextPow2 n
  -- `n .&. (n - 1)` clears the lowest set bit
  | n .&. (n - 1) == 0 = n
  | otherwise = 1 `shiftL` (log2 n + 1)

log2 :: Int -> Int
log2 n = finiteBitSize (0 :: Int) - 1 - countLeadingZeros n


intToRealFloat :: (RealFloat a) => Int -> a
{-# inline intToRealFloat #-}
intToRealFloat = fromIntegral

intToComplex :: (RealFloat a) => Int -> Complex a
{-# inline intToComplex #-}
intToComplex = fromIntegral


{-# inline cmap #-}
cmap :: (Floating a, Unbox a) => (Complex a -> Complex a) -> V.Vector (Complex a) -> V.Vector (Complex a)
cmap = V.map
