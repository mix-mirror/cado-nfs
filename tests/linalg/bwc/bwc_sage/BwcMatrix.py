import os
import re
import sys
import math
import copy
import functools

from sage.matrix.constructor import matrix
from sage.rings.finite_rings.finite_field_constructor import GF
from sage.matrix.special import block_matrix, zero_matrix, identity_matrix
from sage.modules.free_module_element import vector

from .tools import OK, NOK, EXCL
from .tools import u32, s32
from .BwcParameters import BwcParameters
from .BwcBalancing import BwcBalancing, BwcShuffling


class DecorrelatedMatrix(object):
    def __init__(self, parent):
        self.parent = parent
        self.params = parent.params

    def operate(self, y):
        """
        multiplication operator for the decorrelated matrix
        (matrix times vector)
        """

        # TODO: This really tells us that we should not be using
        # operator overload blindly as we do here.
        if isinstance(y, DecorrelatedMatrix):
            raise ValueError("please do not do things like MQ*MQ")

        if self.params.is_nullspace_right():
            return self.parent.operate(self.parent.Q * y)
        else:
            return self.parent.Q.transpose() * self.parent.operate(y)

    def operate_transpose(self, x):
        """
        multiplication operator for the decorrelated matrix
        (vector time matrixvector)
        """

        # TODO: This really tells us that we should not be using
        # operator overload blindly as we do here.
        if isinstance(x, DecorrelatedMatrix):
            raise ValueError("please do not do things like MQ*MQ")

        if self.params.is_nullspace_right():
            return self.parent.Q.transpose() * self.parent.operate_transpose(x)
        else:
            return self.parent.operate_transpose(self.parent.Q * x)

    def __pow__(self, exponent):
        return BwcMatrixPower(self, exponent)


class BwcMatrix(object):
    """
    This class contains stuff to read a matrix in serialized format, but
    also to fetch information related to a given balancing of a matrix,
    or to provide an abstraction of a chain of matrices.

    Sanity checks that compare the on-disk representation with the
    original matrix are included.

    The left and right multiplication operators do not include the
    decorrelating permutation. This means that multplying by this matrix
    is really the same as multiplying by the original matrix that is
    pointed to by the filename, and eventually the output of the whole
    calculation can be checked against this algorithm.

    For cases where the decorrelated matrix M*Q is needed instead, a thin
    wrapper class is provided.

    """
    def __init__(self,
                 params: BwcParameters,
                 matrix=None,
                 wdir=None,
                 balancing_filename=None,
                 multi_matrix=False
                 ):
        if wdir is None:
            self.wdir = os.path.dirname(matrix)
        else:
            self.wdir = wdir
        self.params = params
        self.multi_matrix = multi_matrix

        self.chain = []
        if not multi_matrix:
            self.chain = [BwcOneMatrix(params,
                                       matrix,
                                       wdir,
                                       balancing_filename=balancing_filename)]
        else:
            mm = matrix.split(",")
            if balancing_filename:
                bb = balancing_filename.split(",")
                if len(mm) != len(bb):
                    these = "matrix= and balancing_filename="
                    suck = "should have the same number of items"
                    raise ValueError(f"{these} {suck}")
            else:
                bb = [None for m in mm]
            self.chain = [BwcOneMatrix(params,
                                       mm[i],
                                       wdir,
                                       balancing_filename=bb[i])
                          for i in range(len(mm))]

        self.nrows = None
        self.ncols = None

    @property
    def ncols_orig(self):
        """
        This property is exposed just for convenience, for compatibility
        with the single-matrix case. Note that M.read() must have been
        called first
        """
        return self.chain[-1].ncols_orig

    @property
    def nrows_orig(self):
        """
        This property is exposed just for convenience, for compatibility
        with the single-matrix case. Note that M.read() must have been
        called first
        """
        return self.chain[0].nrows_orig

    def dimensions(self):
        # TODO we should probably get rid of this. There should be no
        # such thing as the dimensions being dependent on
        # self.params.nullspace, really
        assert self.chain[0].nrows is not None
        if self.params.is_nullspace_right():
            return (self.chain[0].nrows, self.chain[-1].ncols)
        else:
            return (self.chain[-1].ncols, self.chain[0].nrows)

    def read(self, force_square=False):

        for c in self.chain:
            c.read()

        if force_square:
            nro = self.chain[0].nrows_orig
            nco = self.chain[-1].ncols_orig
            self.chain[0].pad(outer=(nco, 0), only_extend=False)
            self.chain[-1].pad(outer=(0, nro), only_extend=False)

    def decorrelated(self):
        return DecorrelatedMatrix(self)

    def operate(self, y):
        """
        multiplication operator for this matrix (matrix times vector).
        This multiplication does not include any decorrelated
        permutation.

        y is expected to be a vertical vector with as many rows as the
        dimension of interest of the first matrix that participates in
        the evaluation (i.e., self.ncols, which is the number of columns
        of self.chain[-1] if nullspace=right, and self.nrows, which is
        the number of rows of self.chain[0], otherwise).
        """
        if self.params.is_nullspace_right():
            # we want a right fold, but functools.reduce is a left fold.
            # We can get by with a few reversals, including a reversal of
            # the actual multiplication operation...
            return functools.reduce(lambda v, m: m.operate(v),
                                    reversed(self.chain),
                                    y)
        else:
            # same as: prod(self.chain, y.transpose()).transpose()
            return functools.reduce(lambda v, m: m.operate(v),
                                    self.chain,
                                    y)

    def operate_transpose(self, x):
        """
        multiplication operator for this matrix (vector times matrix),
        This multiplication does not include any decorrelated
        permutation.

        x is expected to be a vertical vector with as many rows as the
        dimension of interest of the first matrix that participates in
        the evaluation. Because we're implementing the transpose operator
        here, this means the number of rows of self.chain[0] if
        nullspace=right, and the number of columns of self.chain[-1]
        otherwise)
        """
        if self.params.is_nullspace_right():
            # same as: prod(self.chain, x.transpose()).transpose()
            return functools.reduce(lambda v, m: m.operate_transpose(v),
                                    self.chain,
                                    x)
        else:
            return functools.reduce(lambda v, m: m.operate_transpose(v),
                                    reversed(self.chain),
                                    x)

    def fetch_balancing(self, nh, nv):
        """
        Based on the filename of the matrix, try to see if a balancing is
        defined, read it, and read all submatrices. Of course this
        assumes that the submatrices were savec during the balancing
        operation done by the C++ code
        """
        for i, c in enumerate(self.chain):
            if (i & 1) == 0:
                c.fetch_balancing(nh, nv)
            else:
                c.fetch_balancing(nv, nh)

        self.Q = self.chain[-1].S.shuf.matrix()

    def check_balancing_submatrices_are_well_formed(self):
        """
        This checks that given the balancing permutations, the
        concatenation of the different submatrices at least gives a
        matrix that is consistent with the number of rows and columns of
        the matrix we supposedly started with.

        It is important to understand that this function does not check
        data that is read by self.read(), and it's on purpose.
        """

        for c in self.chain:
            c.check_balancing_submatrices_are_well_formed()

    def check_balancing_submatrices_consistency_with_M(self):
        """
        This does what the previous code does not do: check that the
        concatenation of the submatrices is consistent with the matrix we
        started with. Of course, to do so, we must have called
        self.read() first.
        """

        for c in self.chain:
            c.check_balancing_submatrices_consistency_with_M()

    def check(self):
        what = "that submatrices are consistent with the balancing"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_are_well_formed()
        print(f"Checking {what} ...  {OK}")

        what = "that submatrices are consistent with the matrix M"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_consistency_with_M()
        print(f"Checking {what} ... {OK}")

    def __pow__(self, exponent):
        return BwcMatrixPower(self, exponent)

    def __str__(self):
        files = ", ".join([x.filename for x in self.chain])
        if self.nrows is not None:
            dims = f"{self.nrows}x{self.ncols} matrix"
            me = f"{dims} from {files}"
        else:
            me = f"matrix from {files} (not read yet)"
        if self.params.is_nullspace_left():
            me += ", operating as vector times matrix"
        else:
            me += ", operating as matrix times vector"
        return me

    def __repr__(self):
        me = [repr(self.params)]
        files = ",".join([x.filename for x in self.chain])
        me.append(f"'{files}'")
        me.append(f"wdir='{self.wdir}'")
        if self.chain[0].balancing is not None:
            bb = ",".join([x.balancing.filename for x in self.chain])
            me.append(f"balancing_filename='{bb}'")
        if len(self.chain) > 1:
            me.append('multi_matrix=True')
        me = ", ".join(me)
        return f"BwcMatrix({me})"



class BwcOneMatrix(object):
    """
    This class contains stuff to read a matrix in serialized format, but
    also to fetch information related to a given balancing of a matrix.

    This type
    Sanity checks that compare the on-disk representation with the
    original matrix are included.

    The left and right multiplication operators do not include the
    decorrelating permutation. This means that multplying by this matrix
    is really the same as multiplying by the original matrix that is
    pointed to by the filename, and eventually the output of the whole
    calculation can be checked against this algorithm.

    For cases where the decorrelated matrix M*Q is needed instead, a thin
    wrapper class is provided.

    """
    def __init__(self,
                 params: BwcParameters,
                 matrix=None,
                 wdir=None,
                 balancing_filename=None,
                 chain_position=0
                 ):
        self.params = params
        self.filename = matrix
        self.chain_position = chain_position
        if wdir is None:
            self.wdir = os.path.dirname(self.filename)
        else:
            self.wdir = wdir

        if balancing_filename is not None:
            self.balancing = BwcBalancing(params, balancing_filename)
        else:
            self.balancing = None

        self.__clear_fields_for_read()
        self.__clear_fields_for_fetch_balancing()

        # The submatrices may be stored in transposed order. This is done
        # because it is convenient for the matmul layers to read them so.
        # Sometimes. In truth, we should perhaps get rid of this
        # complication, I don't really know.
        #
        # Anyway, we have a bit of a problem here. Whether the matrices
        # are stored transposed or not depends on the implementation that
        # got used (because the MM layer "suggests" the good dispatching
        # direction to the layer above). However we don't know what layer
        # is used in this particular code.
        #
        # As a matter of fact, all layers currently in use seem to have
        # MM_DIR0_PREFERS_TRANSP_MULT == 1, which means that if dir==0
        # (nullspace=left), then yes, the layer prefers to see the
        # transposed matrix in ram. And if dir==1, then because 1 xor 1 =
        # 0, then no, it doesn't need that.
        if self.params.is_nullspace_left() and self.params.p == 2:
            self.submatrices_are_transposed = True
        elif self.params.is_nullspace_right() and self.params.p != 2:
            self.submatrices_are_transposed = False
        else:
            raise NotImplementedError("Does the MM layer in this case store matrices in transposed order or not ?")

    def dimensions(self):
        """
        This works only if the matrix has been read already. Returns a
        pair (d1,d2), where d1 is the number of coordinates of vectors X
        that go in a multiplication X*M, and d2 for vectors Y in M*Y.
        Therefore if nullspace=right this is (nrows, ncols) of the
        underlying matrix, and if nullspace=left it's (ncols,nrows).
        """
        assert self.nrows is not None
        if self.params.is_nullspace_right():
            return (self.nrows, self.ncols)
        else:
            return (self.ncols, self.nrows)

    def __clear_fields_for_read(self):
        self.M = None
        self.row_weights = None
        self.row_weights_filename = None
        self.nrows = None
        self.col_weights = None
        self.col_weights_filename = None
        self.ncols = None
        self.ncoeffs = 0
        # (nrows, ncols) are in most cases equal, they correspond to the
        # dimension of a square matrix. They're both equal to
        # max(nrows_orig, ncols_orig), which are the dimensions of the
        # matrix as it is found in the data files on disk.
        self.nrows_orig = None
        self.ncols_orig = None

    def __clear_fields_for_fetch_balancing(self):
        self.balancing = None
        self.S = None
        self.submatrices = None
        self.Mx = None
        self.sigma = None
        self.tau = None
        self.Q = None
        self.xQ = None
        self.P = None
        self.Mt = None

    def __repr__(self):
        ma = f"matrix#{self.chain_position}"
        if self.nrows is not None:
            dims = f"{self.nrows}x{self.ncols} {ma} ({self.ncoeffs} coeffs)"
            return f"{dims} from {self.filename}"
        else:
            return f"{ma} from {self.filename} (not read yet)"

    # Not absolutely sure we want to keep these two operators, in fact
    def operate(self, y):
        """
        this implements the multiplication operator for this matrix. What
        this actually does depends on whether we are going to use this
        matrix for a left or a right operation. For this reason, we
        prefer to avoid operator overloading, as much as possible.

        y is expected to be a vertical vector with as many rows as the
        dimension of interest of the matrix (i.e., its number of columns
        if nullspace=right, its number of rows otherwise).
        """
        if self.params.is_nullspace_right():
            return self.M * y
        else:
            return (y.transpose() * self.M).transpose()

    def operate_transpose(self, x):
        """
        This does the transpose operation of self.operate()
        """
        if self.params.is_nullspace_right():
            return (x.transpose() * self.M).transpose()
        else:
            return self.M * x

    def read(self):
        try:
            self.__read()
        except Exception as e:
            # We're really in bad shape if an exception occurs here.
            # We're not even trying to salvage the BwcMatrix object, as
            # the error is most probably obvious.
            print(f"Exception while reading {self.filename} {NOK}",
                  file=sys.stderr)
            raise e

    def __read(self):
        """
        The parameters nrows and ncols are not meant for general use. We
        only use them when we know the balancing above, and the matrices
        that we're reading are just chunks of the bigger matrix
        """

        print(f"Reading {self.filename}")

        self.__clear_fields_for_read()

        self.row_weights = []
        try:
            fn = re.sub(r"\.bin$", ".rw.bin", self.filename)
            self.nrows_orig = os.stat(fn).st_size // 4
            self.row_weights_filename = fn
            with open(self.row_weights_filename, 'rb') as fm:
                while (w := u32(fm, may_fail=True)) is not None:
                    self.row_weights.append(w)
            assert len(self.row_weights) == self.nrows_orig
        except FileNotFoundError:
            self.row_weights_filename = None

        self.col_weights = []
        try:
            fn = re.sub(r"\.bin$", ".cw.bin", self.filename)
            self.ncols_orig = os.stat(fn).st_size // 4
            self.col_weights_filename = fn
            with open(self.col_weights_filename, 'rb') as fm:
                while (w := u32(fm, may_fail=True)) is not None:
                    self.col_weights.append(w)
            assert len(self.col_weights) == self.ncols_orig
        except FileNotFoundError:
            self.col_weights_filename = None

        inline_data = []
        inline_col_weights = []
        inline_row_weights = []
        if self.params.p == 2:
            with open(self.filename, "rb") as fm:
                i = 0
                j = 0
                try:
                    while (length := u32(fm, may_fail=True)) is not None:
                        for jj in range(length):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                pad = j + 1 - len(inline_col_weights)
                                inline_col_weights += [0] * pad
                            inline_col_weights[j] += 1
                            inline_data.append((i, j))
                        i += 1
                        inline_row_weights.append(length)
                except IndexError:
                    what = f"entry {i},{j} in matrix"
                    where = f"at byte 0x{fm.tell():04x} in matrix file"
                    raise ValueError(f"Cannot set {what} {where} {NOK}")
        else:
            with open(self.filename, "rb") as fm:
                i = 0
                j = 0
                v = 0
                try:
                    while (length := u32(fm, may_fail=True)) is not None:
                        for jj in range(length):
                            j = u32(fm)
                            if j >= len(inline_col_weights):
                                pad = j + 1 - len(inline_col_weights)
                                inline_col_weights += [0] * pad
                            v = s32(fm)
                            inline_col_weights[j] += 1
                            inline_data.append((i, j, v))
                        i += 1
                        inline_row_weights.append(length)
                except IndexError:
                    what = f"entry {i},{j} in matrix"
                    where = f"at byte 0x{fm.tell():04x} in matrix file"
                    raise ValueError(f"Cannot set {what} {where} {NOK}")

        if self.row_weights:
            assert self.row_weights == inline_row_weights
        else:
            self.row_weights = inline_row_weights
            self.nrows_orig = len(self.row_weights)

        if self.col_weights:
            icw = inline_col_weights
            assert self.col_weights[:len(icw)] == icw
            assert vector(self.col_weights[len(icw):]) == 0
        else:
            self.col_weights = inline_col_weights
            self.ncols_orig = len(self.col_weights)

        self.ncoeffs = len(inline_data)
        r2 = sum([float(x*x) for x in inline_row_weights]) / self.nrows_orig
        c2 = sum([float(x*x) for x in inline_col_weights]) / self.ncols_orig
        rmean = float(self.ncoeffs / self.nrows_orig)
        rsdev = math.sqrt(r2-rmean**2)
        cmean = float(self.ncoeffs / self.ncols_orig)
        csdev = math.sqrt(c2-cmean**2)

        rowscols = f"{self.nrows_orig} rows {self.ncols_orig} cols"
        coeffs = f"{self.ncoeffs} coefficients"
        stats = f"row: {rmean:.2f}~{rsdev:.2f}, col: {cmean:.2f}~{csdev:.2f}"
        print(f"{rowscols}, {coeffs} ({stats})")

        # prior to any padding, we have a matrix with dimensions
        # nrows_orig and ncols_orig

        self.M = matrix(GF(self.params.p),
                        self.nrows_orig,
                        self.ncols_orig,
                        sparse=True)
        if self.params.p == 2:
            for i, j in inline_data:
                if self.M[i, j] != 0:
                    print(f"Warning: repeated entry in matrix data {EXCL}")
                self.M[i, j] += 1
        else:
            for i, j, v in inline_data:
                if self.M[i, j] != 0:
                    print(f"Warning: repeated entry in matrix data {EXCL}")
                self.M[i, j] += v

        assert self.balancing is None

        if self.balancing is not None:
            self.balancing.read()
            self.S = BwcShuffling(self.params, self)

    def pad(self, outer=None, only_extend=True, force_square=False):
        self.nrows = self.nrows_orig

        # Attention. It may be that the matrix has some zero cols at the
        # end, so if self.ncols_orig was just guessed based on the max
        # column index, we definitely need some correct padding.
        self.ncols = self.ncols_orig

        if outer is not None:
            nro, nco = outer
            if only_extend:
                assert self.nrows_orig <= nro or nro == 0
                assert self.ncols_orig <= nco or nco == 0
            self.nrows = max(nro, self.nrows)
            self.ncols = max(nco, self.ncols)

        if force_square:
            print("Padding to a square matrix")
            self.nrows = max(self.nrows, self.ncols)
            self.ncols = max(self.nrows, self.ncols)

        newrows = self.nrows - self.nrows_orig
        newcols = self.ncols - self.ncols_orig

        if newrows:
            who = f"matrix {self.filename}"
            what = f"{newrows} new rows"
            reach = f"so that we reach {self.nrows} in total"
            print(f"Padding {who} with {what} ({reach})")

        if newcols:
            who = f"matrix {self.filename}"
            what = f"{newcols} new cols"
            reach = f"so that we reach {self.ncols} in total"
            print(f"Padding {who} with {what} ({reach})")

        self.M = block_matrix(2,2,[
            self.M, matrix(self.nrows_orig, newcols),
            matrix(newrows, self.ncols_orig), matrix(newrows, newcols)])
        self.M.subdivide()


    def __subM(self, i, j):
        return self.submatrices[i][j].M

    def fetch_balancing(self, nh, nv):
        """
        Based on the filename of the matrix, try to see if a balancing is
        defined, read it, and read all submatrices. Of course this
        assumes that the submatrices were savec during the balancing
        operation done by the C++ code
        """

        self.__clear_fields_for_fetch_balancing()

        dir, base = os.path.split(self.filename)
        dir = self.wdir
        sdir = base.removesuffix(".bin") + f".{nh}x{nv}"
        bfile = os.path.join(dir, sdir, sdir + ".bin")
        if not os.path.exists(bfile):
            return FileNotFoundError(bfile)
        self.balancing = BwcBalancing(self.params, bfile)
        print("Reading balancing from " + bfile)
        self.balancing.read()
        self.S = BwcShuffling(self.params, self)
        self.submatrices = [[None for j in range(nv)] for i in range(nh)]
        bal = self.balancing
        for i in range(nh):
            for j in range(nv):
                # id = f"{self.balancing.checksum:08x}.h{i}.v{j}"
                id = f"h{i}.v{j}"
                fname = os.path.join(dir, sdir, f"{sdir}.{id}.bin")
                nrp = bal.tr // bal.nh
                ncp = bal.tc // bal.nv
                self.submatrices[i][j] = BwcOneMatrix(self.params, fname)
                if self.submatrices_are_transposed:
                    self.submatrices[i][j].read()
                    self.submatrices[i][j].pad(outer=(ncp, nrp))
                else:
                    self.submatrices[i][j].read()
                    self.submatrices[i][j].pad(outer=(nrp, ncp))
        # Now we want to check the consistency of the matrix that we just
        # read.

        t = lambda x: x  # noqa: E731
        if self.submatrices_are_transposed:
            t = lambda x: x.transpose()  # noqa: E731

        self.Mx = block_matrix(nh, nv,
                               [t(self.__subM(i, j))
                                for i in range(nh)
                                for j in range(nv)])

        # sigma is applied on rows.
        # By convention, the action of sigma is transpose(sigma) on the
        # left. (so that we can have sigma==tau for square matrices)
        if self.S.sr is not None:
            self.sigma = self.S.sr.matrix()
        else:
            self.sigma = identity_matrix(GF(self.params.p), self.balancing.tr)

        # tau is applied on columns
        if self.S.sc is not None:
            self.tau = self.S.sc.matrix()
        else:
            self.tau = identity_matrix(GF(self.params.p), self.balancing.tc)

        if self.S.shuf is not None:
            self.Q = self.S.shuf.matrix()
            self.xQ = self.S.xshuf.matrix()
        else:
            self.Q = identity_matrix(GF(self.params.p), self.balancing.nc)
            self.xQ = identity_matrix(GF(self.params.p), self.balancing.tc)

        self.P = self.S.pr.matrix()

    def check_balancing_submatrices_are_well_formed(self):
        """
        This checks that given the balancing permutations, the
        concatenation of the different submatrices at least gives a
        matrix that is consistent with the number of rows and columns of
        the matrix we supposedly started with.

        It is important to understand that this function does not check
        data that is read by self.read(), and it's on purpose.
        """

        bal = self.balancing

        print(f"Checking {self.filename} ({bal.nh}x{bal.nv} balancing)")

        A = (self.P*self.sigma).transpose() * self.Mx * self.tau

        what = "Reconstructed matrix from submatrices"
        if A[bal.nr:bal.tr] != 0:
            raise ValueError(f"{what} has garbage trailing rows {NOK}")
        if A.transpose()[bal.nc:bal.tc] != 0:
            raise ValueError(f"{what} has garbage trailing columns {NOK}")

    def check_balancing_submatrices_consistency_with_M(self):
        """
        This does what the previous code does not do: check that the
        concatenation of the submatrices is consistent with the matrix we
        started with. Of course, to do so, we must have called
        self.read() first.
        """

        bal = self.balancing

        if self.M is None:
            raise ValueError(f".read() must be called first {EXCL}")

        if bal is None:
            raise ValueError(f".fetch_balancing() must be called first {EXCL}")

        if self.nrows_orig != bal.nr or self.ncols_orig != bal.nc:
            what = "inconsistent with number of rows/columns of the matrix"
            raise ValueError(f"{self.balancing.filename} is {what} {NOK}")

        A = (self.P*self.sigma).transpose() * self.Mx * self.tau

        B = self.M * self.Q
        sub = lambda T: T.submatrix(0, 0, self.nrows, self.ncols)  # noqa: E731

        what = "Reconstructed matrix from submatrices"
        if sub(A) != B:
            raise ValueError(f"{what} does not match M " + NOK)

        self.Mt = zero_matrix(bal.tr, bal.tc)
        self.Mt[:self.nrows, :self.ncols] = self.M

        if A != self.Mt * self.xQ:
            # It's not really possible for this to happen if the above
            # check has passed, so let's just give the same error message
            # anyway.
            raise ValueError(f"{what} does not match M " + NOK)

    def check(self):
        nh = self.balancing.nh
        nv = self.balancing.nv

        what = f"that submatrices are well formed for {nh}*{nv} balancing"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_are_well_formed()
        print(f"Checking {what} ...  {OK}")

        what = "that submatrices are consistent with the matrix M"
        print(f"Checking {what} ...")
        self.check_balancing_submatrices_consistency_with_M()
        print(f"Checking {what} ... {OK}")

    def __pow__(self, exponent):
        return BwcMatrixPower(self, exponent)


class BwcMatrixPower(object):
    def __init__(self, M, exponent):
        self.M = M
        self.exponent = exponent

    def operate(self, y):
        v = copy.copy(y)
        for i in range(self.exponent):
            v = self.M.operate(v)
        return v

    def operate_transpose(self, y):
        v = copy.copy(y)
        for i in range(self.exponent):
            v = self.M.operate_transpose(v)
        return v
