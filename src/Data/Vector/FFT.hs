{-# language BangPatterns        #-}
{-# language LambdaCase          #-}

module Data.Vector.FFT where

import Control.Monad (when)
import Control.Monad.Primitive (PrimMonad(..))

import Control.Monad.ST (runST)
import Data.Bits (shiftR,shiftL,(.&.),(.|.))
import Data.Bool (Bool,otherwise)
import Data.Complex (Complex(..),conjugate)

import Data.Vector.Unboxed as V (Vector, Unbox(..), map, length, unsafeFreeze)
import qualified Data.Vector.Unboxed.Mutable as VM (read, write, copy, new)


import Prelude hiding (read)

{-# RULES
"fft/ifft" forall x. fft (ifft x) = x
"ifft/fft" forall x. ifft (fft x) = x
  #-}

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two.

{-# inlinable [1] fft #-}
fft arr = if arrOK arr
  then runST $ do {
      marr <- copyWhole arr
    ; mfft marr
    ; V.unsafeFreeze marr
  }
  else Prelude.error "Data.Vector.FFT.fft: bad array length"

-- | Inverse fast Fourier transform.

{-# inlinable [1] ifft #-}
ifft arr = if arrOK arr
  then
    let lenComplex = intToComplexDouble (V.length arr)
    in cmap ((/ lenComplex) . conjugate) . fft . cmap conjugate $ arr
  else error "Data.Vector.FFT.ifft: bad vector length"


{-# inline copyWhole #-}
copyWhole arr = do
  let len = V.length arr
  marr <- VM.new len
  VM.copy marr arr
  pure marr


{-# inline arrOK #-}
arrOK arr =
  let n = V.length arr
  in (1 `shiftL` log2 n) == n

-- | Radix-2 decimation-in-time fast Fourier Transform.
--   The given array must have a length that is a power of two,
--   though this property is not checked.

mfft mut = do {
    let len = V.length mut
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
                                ; let !c = cos a
                                      !s = sin a
                                      d = (c*xi1 - s*yi1) :+ (s*xi1 + c*yi1)
                                ; ci <- VM.read mut i
                                ; VM.write mut i1 (ci - d)
                                ; VM.write mut i (ci + d)
                                ; butterfly (i + l2)
                              }
                      ; butterfly j
                    }
            ; flight 0 0
         }
  ; bitReverse 0 0
}

-- wildcard cases should never happen. if they do, really bad things will happen.
b,s :: Int -> Int
b = \case { 0 -> 0x02; 1 -> 0x0c; 2 -> 0xf0; 3 -> 0xff00; 4 -> wordToInt 0xffff0000; 5 -> wordToInt 0xffffffff00000000; _ -> 0; }
s = \case { 0 -> 1; 1 -> 2; 2 -> 4; 3 -> 8; 4 -> 16; 5 -> 32; _ -> 0; }
{-# inline b #-}
{-# inline s #-}

log2 :: Int -> Int
log2 v0 = if v0 <= 0
  then error $ "Data.Vector.FFT: nonpositive input, got " ++ show v0
  else go 5 0 v0
  where
    go !i !r !v
      | i == -1 = r
      | v .&. b i /= 0 =
          let si = s i
          in go (i - 1) (r .|. si) (v `shiftR` si)
      | otherwise = go (i - 1) r v


{-# inline swap #-}
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

wordToInt :: Word -> Int
{-# inline wordToInt #-}
wordToInt = fromIntegral

intToComplexDouble :: Int -> Complex Double
{-# inline intToComplexDouble #-}
intToComplexDouble = fromIntegral


{-# inline cmap #-}
cmap = V.map
