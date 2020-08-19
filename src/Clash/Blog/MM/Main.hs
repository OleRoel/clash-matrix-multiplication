
module Clash.Blog.MM.Main where

import Clash.Prelude

type Matrix m n a = Vec m (Vec n a)

emptyMatrix :: a -> Matrix m n a
emptyMatrix = repeat . repeat

nullMatrix :: Num a => Matrix m n a
nullMatrix = repeat (repeat 0)

-- | Dot product
dot
  :: KnownNat n
  -- ^ Constraint 1: Store length information at runtime too
  => 1 <= n
  -- ^ Constraint 2: Vectors must be at least of length one
  => Vec n Int
  -> Vec n Int
  -> Int
dot vec1 vec2 =
  sum (zipWith (*) vec1 vec2)

-- | Matrix/vector multiplication
mvMult
  :: KnownNat n
  => 1 <= n
  -- ^ Constraints needed for `dot`
  => Matrix m n Int
  -- ^ Matrix with `m` rows, `n` columns
  -> Vec n Int
  -- ^ Vector with `n` integers
  -> Vec m Int
mvMult mat vec =
  map (dot vec) mat

-- | Matrix/matrix multiplication
mmMult
  :: an ~ bm
  -- ^ Number of columns of matrix A must be
  -- equal to the number of rows in matrix B.
  => 1 <= bm
  => KnownNat bn
  => KnownNat bm
  => Matrix am an Int
  -> Matrix bm bn Int
  -> Matrix am bn Int
mmMult mat1 mat2 =
  map (mvMult (transpose mat2)) mat1

mmmult2d
  :: forall a_m a_n b_m b_n aa_m aa_n bb_m bb_n aa_sm aa_sn bb_sm bb_sn
  -- ^ Explicit definition of type variables in order to use them in function body

      -- Clock and reset lines for registers
  . ( SystemClockResetEnable

    , KnownNat aa_m
    , KnownNat aa_n
    , KnownNat bb_m
    , KnownNat bb_n
    , KnownNat aa_sm
    , KnownNat aa_sn
    , KnownNat bb_sm
    , KnownNat bb_sn

    -- Enforce proper matrix dimensions:
    , 1 <= a_m
    , 1 <= a_n
    , 1 <= b_m
    , 1 <= b_n
    , a_n ~ b_m

    -- Constrain submatrices:
    , a_m ~ (aa_m * aa_sm)
    , a_n ~ (aa_n * aa_sn)
    , b_m ~ (bb_m * bb_sm)
    , b_n ~ (bb_n * bb_sn)
    , 1 <= aa_sm
    , 1 <= aa_sn
    , 1 <= bb_sm
    , 1 <= bb_sn
    , bb_sm ~ aa_sn )

  -- Allow user to pass submatrix sizes:
  => SNat aa_m
  -- ^ Number of rows in submatrix of AA
  -> SNat aa_sn
  -- ^ Number of columns in submatrix of AA
  -> SNat bb_n
  -- ^ Number of columns in submatrix of BB

  -- Matrices to multiply:
  -> Signal System (Maybe (Matrix a_m a_n Int, Matrix b_m b_n Int))

  -- Result returned after calculating for a while:
  -> Signal System (Maybe (Matrix a_m b_n Int))
mmmult2d aa_sn aa_m bb_n ab =
  mealy mmmult2dmealy state ab'
    where
      -- Take input matrices, and split them into smaller ones. The outer fmap
      -- maps over each value in the signal, the inner fmap applies the function
      -- `splitab` on the inner value of Maybe (if it is not Nothing).
      ab' = fmap (fmap splitab) ab

      -- Initial state for mealy machine:
      state =
        ( Nothing                -- No matrices saved yet
        , minBound               -- Counter at zero
        , emptyMatrix nullMatrix -- Matrix with zero-matrices
        )

      -- Split matrices into matrix with submatrices
      splitab (a, b) =
        ( msplit a :: Matrix aa_m aa_n (Matrix aa_sm aa_sn Int)
        , msplit b :: Matrix bb_m bb_n (Matrix bb_sm bb_sn Int)
        )

-- | Same as (!!) but guaranteed to succeed as any
-- value in `Index n` can never exceed `n-1`.
index
  :: KnownNat n
  => Vec n a
  -> Index n
  -> a
index = (!!)

-- | mmmult2dmealy describes a single caclulation step. It returns a result only
-- when it's ready. To be used as mealy machine.
mmmult2dmealy (Nothing, _, _) Nothing =
  -- No input nor state, do nothing:
  ((Nothing, minBound, emptyMatrix nullMatrix), Nothing)

mmmult2dmealy _ matrices@(Just _) =
  -- Input; reset progress so far (if any)
  ((matrices, minBound, emptyMatrix nullMatrix), Nothing)

mmmult2dmealy (matrices@(Just (matrixAA, matrixBB)), counter, matrixRR) _ = (state', output)
  -- Continue calculating, return result if ready
  where
    -- If we're at the counter's maximum, we're done after this cycle
    done = counter == maxBound

    -- Increase counter tuple by one. Wrap around if maximum is reached.
    counter' = succWrap counter

    -- Calculate new state; if we're done, reset it.
    state'
      | done      = (Nothing,  counter', emptyMatrix nullMatrix)
      | otherwise = (matrices, counter', matrixRR')

    -- Output only if we're done calculating
    output
      | done      = Just (mmerge matrixRR')
      | otherwise = Nothing

    -- Determine order of fetching from A or B and storing it in R.
    (aColI,  _,     aRowI) = counter
    (bRowI,  bColI, _    ) = counter
    (_,      rColI, rRowI) = counter

    -- Fetch submatrices and partial result
    subA = (matrixAA `index` aRowI) `index` aColI
    subB = (matrixBB `index` bRowI) `index` bColI
    subR = (matrixRR `index` rRowI) `index` rColI

    -- Calculate new partial result, store it in matrix R
    subR'     = madd subR (mmMult subA subB)
    matrixRR' = replaceMatrixElement matrixRR rRowI rColI subR'

main :: IO ()
main = do
  let a = (1 :> 2 :> Nil)
  let b = (5 :> 7 :> Nil)
  putStrLn $ show $ dot a b
