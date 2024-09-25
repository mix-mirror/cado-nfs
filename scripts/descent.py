#!/usr/bin/env python3
#
# Examples of parameters can be found in the parameters/dlp directory

# Algorithms
#   UpperClass:
#       This is the initialization of the descent. An extended gcd is
#       done between the target and p to get two "rational
#       reconstructions", and then a sieving procedure is done to find a
#       linear combination of these that is smooth (idea taken from
#       Joux-Lercier). This sieving is done with las, with two linear
#       polynomials.
#       At the end, we have target = num/den, where num and den are
#       smooth, with a smoothness bound that is large than the one that
#       was used in the sieving / linear algebra.
#   MiddleClass:
#       For all the primes dividing num and den that are larger than the
#       large prime bound, a "special-q descent" is performed, in order
#       to rewrite them in terms of smaller and smaller elements, until
#       everything is known. This is done with the las_descent program.
#   LowerClass:
#       This step is just putting everything together.
#       In practice, this means computing appropriate Schirokauer maps
#       and propagating the known logarithms in relations in order to
#       deduce the unknown ones. The main tool for that is the
#       reconstructlog program. Some ugly modifications of input files
#       are necessary.


# TODO: do we want to have default values, as done here, for the --init-*
# arguments ? I would say no.

# TODO: of course this is currently kludgy, and does not properly wield
# the power of the cadofactor python programs. I need to understand that
# stuff better.

# TODO: keep las awake, if needed.

# TODO: make output less verbose.

import os
import io
import sqlite3
import subprocess
import sys
import argparse
import re
import math
import time
import tempfile
import shutil
import functools
import itertools
import random
if sys.version_info >= (3,8):
    from descent_helper_asyncio import monitor_important_files
else:
    from descent_helper_fallback import monitor_important_files


has_hwloc = None

# This gives the boring info about the file names for everything, and
# the boilerplate arguments to be fed to binaries.
class GeneralClass(object):

    def declare_args(parser):
        # parser.add_argument("--no-wipe",
        #         help="Keep working files",
        #         action="store_true")
        parser.add_argument("--sm-mode",
                help="Select SM mode",
                type=str)
        parser.add_argument("--datadir",
                help="cadofactor working directory",
                type=str)
        parser.add_argument("--prefix",
                help="project prefix",
                type=str)
        parser.add_argument("--db",
                help="SQLite db file name",
                type=str)
        # the next few are optional file names
        parser.add_argument("--tmpdir", help="Temporary working directory")
        parser.add_argument("--cadobindir",
                help="Cado build directory",
                required=True,
                type=str)
#        parser.add_argument("--todofile",
#                help="Output primes to this toto file",
#                type=str)
        parser.add_argument("--poly",
                help="Polynomial file",
                type=str)
        parser.add_argument("--renumber",
                help="Renumber file",
                type=str)
        parser.add_argument("--fb1",
                help="Factor base file for the algebraic side",
                type=str)
        parser.add_argument("--log",
                help="File with known logs",
                type=str)
        parser.add_argument("--gfpext",
                help="Degree of extension (default 1)",
                type=int)
        # This one applies to both las in the initial step
        parser.add_argument("--threads",
                help="Number of threads to use",
                type=int, default=4)
        # the arguments below are really better fetched from the
        # database.
        parser.add_argument("--ell", help="Group order (a.k.a. ell)")
        parser.add_argument("--nsm0", help="Number of Schirokauer maps on side 0")
        parser.add_argument("--nsm1", help="Number of Schirokauer maps on side 1")
        # Those are used both for the middle and lower levels of the
        # descent.
        for side in range(2):
            parser.add_argument("--lpb%d" % side,
                    help="Default large prime bound on side %d" % side,
                    required=True,
                    type=int)

    def __init__(self, args):
        self._conn = None
        self.args = args
        if bool(args.db) == bool(args.prefix and args.datadir):
            raise ValueError("Either --db (with an sqlite db) or the combo --prefix + --datadir must be specified")
        if args.tmpdir:
            self._tmpdir = args.tmpdir
            # do mkdir ???
        else:
            self._tmpdir = tempfile.mkdtemp(dir="/tmp")
        self.hello()
        self.__load_badidealdata()
        self.logDB = LogBase(self)
        self.initrandomizer = 1
    
    def __connect(self):
        if args.db and not self._conn:
            self._conn = sqlite3.connect(args.db)

    def __getdb(self, query):
        if not args.db:
            return None
        self.__connect()
        self._cursor = self._conn.cursor()
        self._cursor.execute(query)
        v=self._cursor.fetchone()
        self._cursor.close()
        return v

    def __getfile(self, shortname, typical, table, key):
        try:
            v=self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if args.db and table:
            v=self.__getdb("select value from %s where kkey='%s'" % (table, key))
            if v is not None and len(v) > 0:
                return os.path.join(os.path.dirname(args.db), v[0])
        elif args.datadir and args.prefix:
            return os.path.join(args.datadir, args.prefix + "." + typical)
        raise ValueError("no %s file known" % shortname)

    def __getarg(self, shortname, table, key):
        try:
            v=self.args.__dict__[shortname]
            if v:
                return v
        except KeyError:
            pass
        if args.db:
            v=self.__getdb("select value from %s where kkey='%s'" % (table, key))
            if v is not None and len(v) > 0:
                return v[0]
        raise ValueError("no %s parameter known" % shortname)

    def prefix(self):
        if args.prefix:
            return args.prefix
        else:
            return os.path.basename(args.db).split('.')[0]

    def datadir(self):
        if args.datadir:
            return args.datadir
        elif args.db:
            return os.path.dirname(args.db)
        else:
            raise ValueError("Need --datadir or --db with an sqlite db")

    def poly(self):
        return self.__getfile("poly", "poly", "polyselect2", "polyfilename")
    def renumber(self):
        return self.__getfile("renumber", "renumber.gz", "freerel", "renumberfilename")

    def log(self):
        return self.__getfile("log", "dlog", "reconstructlog", "dlog")
    def badideals(self):
        return os.path.join(self.datadir(), self.prefix() + ".badideals")
    def badidealinfo(self):
        return os.path.join(self.datadir(), self.prefix() + ".badidealinfo")
    def fb1(self):
        return self.__getfile("fb1", "roots1.gz", "factorbase", "outputfile")
    def fb0(self):
        return self.__getfile("fb0", "roots0.gz", "factorbase", "outputfile")
    def ell(self):
        return int(args.ell)
    def lpb0(self):
        return args.lpb0
    def lpb1(self):
        return args.lpb1
    def tmpdir(self):
        return self._tmpdir
    def threads(self):
        return int(args.threads)
    def poly_data(self):
        d={}
        with open(self.poly(), "r") as file:
            for line in file:
                if re.match("^\s*#", line):
                    continue
                if re.match("^\s*$", line):
                    continue
                key,value=line.split(":")
                key = key.strip()
                foo = re.match("^([cY])(\d+)$", key)
                if foo:
                    s,i=foo.groups()
                    if s not in d:
                        d[s]=[]
                    while int(i) >= len(d[s]):
                        d[s]+=[None]
                    d[s][int(i)]=value.strip()
                else:
                    d[key] = value.strip()
            if 'poly0' in d:
                assert 'Y' not in d
                d['Y'] = [ int(x) for x in d["poly0"].split(',') ]
            if 'poly1' in d:
                assert 'c' not in d
                d['c'] = [ int(x) for x in d["poly1"].split(',') ]
        return d

    def p(self):
        d=self.poly_data()
        return int(d["n"])
    def extdeg(self):
        if args.gfpext:
            return args.gfpext
        else:
            return 1

    def target(self):
        if self.extdeg() == 1:
            return int(args.target)
        else:
            return [int(x) for x in args.target.split(",")]
    
    # short name for the target, to be used in filenames
    def short_target(self):
        target = str(self.args.target)
        if len(target) <= 20:
            return target
        else:
            return target[:10] + "..." + target[-10:]

    def has_rational_side(self):
        d=self.poly_data()
        return len(d["Y"]) == 2
    def rational_poly():
        d=self.poly_data()
        assert len(d["Y"]) == 2
        return [int(x) for x in d["Y"]]
    def algebraic_poly():
        d=self.poly_data()
        return [int(x) for x in d["c"]]

    def cleanup(self):
        if False:
#        if not self.args.tmpdir and not self.args.no_wipe:
            shutil.rmtree(self.tmpdir())

    def __del__(self):
        if self._conn:
            self._conn.close()

    def descentinit_bin(self):
        return os.path.join(args.cadobindir, "misc", "descent_init_Fp")

    def las_bin(self):
        return os.path.join(args.cadobindir, "sieve", "las")
    def sm_simple_bin(self):
        return os.path.join(args.cadobindir, "filter", "sm_simple")
    def numbertheory_bin(self):
        return os.path.join(args.cadobindir, "utils", "numbertheory_tool")

    def lasMiddle_base_args(self):
        # TODO add threads once it's fixed.
        s=[
            self.las_bin() + "_descent",
            "--recursive-descent",
            "--allow-largesq",
            "--never-discard",  # useful for small computations.
            "--renumber", self.renumber(),
            "--log", self.log(),
            "--fb1", self.fb1(),
            "--poly", self.poly(),
          ]
        if not self.has_rational_side():
            s.append("--fb0")
            s.append(self.fb0())
        return [ str(x) for x in s ]

    # There's no las_init_base_args, since DescentUpperClass uses only
    # its very own arguments.

    def hello(self):
        print("Working in GF(p), p=%d" % self.p())
        print("Subgroup considered in GF(p)^* has size %d" % self.ell())
        print("prefix is %s" % self.prefix())
        errors=[]
        if not os.path.exists(self.las_bin()):
            errors.append("las not found (make las ?)")
        if not os.path.exists(self.las_bin() + "_descent"):
            errors.append("las_descent not found (make las_descent ?)")
        if not os.path.exists(self.sm_simple_bin()):
            errors.append("sm_simple not found (make sm_simple ?)")
        if not os.path.exists(self.numbertheory_bin()):
            errors.append("numbertheory_tool not found (make numbertheory_tool ?)")
        for f in [ self.log(), self.poly(), self.renumber(), self.log(), self.fb1() ]:
            if not os.path.exists(f):
                errors.append("%s missing" % f)
        if len(errors):
            msg = "Some data files and/or binaries missing:\n"
            msg += "\n".join(["\t"+x for x in errors])
            raise RuntimeError(msg)

    # self.list_badideals will contain a list of (p,r,side)
    # self.list_bad_ncols will contain a list of the corresponding nb of cols
    # self.badidealdata will contain a list of
    #     (p, k, rk, side, [exp1, exp2, ..., expi])
    def __load_badidealdata(self):
        # Note that we can as well get this information from the renumber
        # file.
        if not os.path.exists(self.badideals()) or not os.path.exists(self.badidealinfo()):
            call_that = [ self.numbertheory_bin(),
                            "-poly", self.poly(),
                            "-badideals", self.badideals(),
                            "-badidealinfo", self.badidealinfo(),
                            "-ell", self.ell()
                        ]
            call_that = [str(x) for x in call_that]
            print("command line:\n" + " ".join(call_that))
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(call_that, stderr=devnull)

        self.list_badideals = []
        self.list_ncols = []
        with open(self.badideals(), 'r') as bad:
            for line in bad:
                if line[0] == '#':
                    continue
                foo = re.match("^(\d+),(\d+):(\d+): (\d+)$", line)
                if foo:
                    self.list_badideals.append((int(foo.groups()[0]),
                        int(foo.groups()[1]), int(foo.groups()[2])))
                    self.list_ncols.append(int(foo.groups()[3]))
                else:
                    raise ValueError("Error while reading %s" % self.badideal())

        self.badidealdata = []
        with open(self.badidealinfo(), 'r') as bad:
            for line in bad:
                if line[0] == '#':
                    continue
                pattern = "^(\d+) (\d+) (\d+) (\d+) (.+)$"
                foo = re.match(pattern, line)
                if foo:
                    self.badidealdata.append((
                        int(foo.groups()[0]), # p
                        int(foo.groups()[1]), # k
                        int(foo.groups()[2]), # rk
                        int(foo.groups()[3]), # side
                        [ int(x) for x in foo.groups()[4].split() ] # exp
                        ))
                else:
                    raise ValueError("Error while reading %s" %
                        self.badidealinfo())
        print ("Bad ideal information loaded: %s bad ideals, and %s lines in badidealinfo"
                % (str(len(self.list_badideals)), str(len(self.badidealdata))))
        print ("badideal data: %s" % str(self.badidealdata))


def check_result(two, log2, z, logz, p, ell):
    assert (p-1) % ell == 0
    assert pow(z, log2*((p-1) // ell), p) == pow(2, logz*((p-1) // ell), p)
    print ("Final consistency check ok!")


# A memory image of the reconstructlog.dlog file.
class LogBase(object):
    def __init__(self, general):
        self.known={}
        self.badideals=[]
        self.SMs = [ [], [] ]
        self.fullcolumn = None
        try:
            print ("--- Reading %s to find which are the known logs ---" % general.log())
            def process(line):
                index,p,side,r,*value = line.split()
                if p == "bad" and side == "ideals":
                    self.badideals.append(int(r))
                elif p == "SM":
                    self.SMs[int(side)].append(int(value[0]))
                else:
                    # for rational side, we actually don't have the root.
                    # We force it to be on side 0 (legacy, again...)
                    if r == "rat":
                        assert int(side) == 0;
                        r = -1
                    else:
                        r = int(r, 16)
                    self.known[(int(p, 16),r,int(side))] = int(value[0])

            with open(general.log(),'r') as file:
                line = file.readline()
                m = re.match("^(\w+) added column (\d+)$",line)
                if m:
                    self.fullcolumn = int(m.groups()[1])
                else:
                    self.fullcolumn = None
                    process(line)
                for i,line in enumerate(file):
                    if i % 1000000 == 0:
                        print("Reading line %d" % i)
                    process(line)
            print("Found %d bad ideals, %d known logs, and %d,%d SMs in %s"
                %(len(self.badideals), len(self.known),
                len(self.SMs[0]), len(self.SMs[1]), general.log()))
        except:
            raise ValueError("Error while reading %s" % general.log())
    def has(self,p,r,side):
        return (p,r,side) in self.known
    def get_log(self, p,r,side):
        if general.has_rational_side() and side == 0:
            r = -1;
        if (p,r,side) in self.known:
            return self.known[(p,r,side)]
        else:
            return None
    def add_log(self,p,r,side,log):
        if general.has_rational_side() and side == 0:
            r = -1;
        self.known[(p,r,side)] = log
    def bad_ideal(self,i):
        return self.badideals[i]
    def allSM(self,i):
        SM = self.SMs[0] + self.SMs[1]
        return SM[i]
    def nSM(self,side):
        return len(self.SMs[side])
    def SM(self,side,i):
        return self.SMs[side][i]
    def full_column(self):
        return self.fullcolumn

def a_over_b_mod_p(a, b, p):
    if b%p == 0:
        return p
    ib = pow(b, p-2, p)
    return (a*ib) % p

def is_a_over_b_equal_r_mod_pk(a, b, rk, p, pk):
    if b%p != 0:   # non-projective
        if rk >= pk:
            return False
        return (a-b*rk) % pk == 0
    else: # projective
        if rk < pk:
            return False
        return (b-a*rk) % pk == 0

class ideals_above_p(object):
    def __is_badideal(self, p, r, side, general):
        return (p,r,side) in general.list_badideals

    def __handle_badideal(self, p, k, a, b, r, side, general):
        baddata = general.badidealdata
        expo = []
        badid = None
        for X in baddata:
            if side != X[3] or p != X[0]:
                continue;
            pk = pow(p, X[1])
            rk = X[2]
            if not is_a_over_b_equal_r_mod_pk(a,b,rk,p,pk):
                continue;
            vals = X[4]
            badid = (p, r, side)
            for v in vals:
                if v > 0:
                    exp = v
                else:
                    assert k >= -v
                    exp = k+v
                expo.append(exp)
            break
        if badid == None:
            raise ValueError("Error while handling badideal p=%d side=%d a=%d b=%d" % (p, side, a, b))
        return {"badid":badid, "exp":expo}

    def __init__(self, p, k, a, b, side, general):
        self.logDB = general.logDB
        self.p = p
        self.k = k
        self.side = side
        if general.has_rational_side() and side == 0:
            self.r = -1
        else:
            self.r = a_over_b_mod_p(a, b, p)
        self.isbad = self.__is_badideal(p, self.r, side, general)
        if self.isbad:
            self.bads = self.__handle_badideal(p, k, a, b, self.r, side, general)

    # return an unreduced virtual log or None if unknown.
    def get_log(self):
        if not self.isbad:
            l = self.logDB.get_log(self.p, self.r, self.side)
            if l == None:
                return None
            else:
                return self.k * l
        else:
            ind_b = 0
            ind_l = 0
            while general.list_badideals[ind_b] != (self.p, self.r, self.side):
                ind_l += general.list_ncols[ind_b]
                ind_b += 1
            logs = [self.logDB.bad_ideal(x) for x in range(ind_l,
                ind_l+general.list_ncols[ind_b]) ]
            assert len(logs) == len(self.bads["exp"])
            log = 0
            for i in range(len(logs)):
                log = log + logs[i]*self.bads["exp"][i]
            return log

class object_holder():
    def __init__(self, v=None):
        self.v = v
    def set(self, v):
        self.v = v
    def unset(self):
        self.v = None

class DescentUpperClass(object):
    def declare_args(parser):
        c = " (specific for the descent bootstrap)"
        parser.add_argument("--init-tkewness",
                help="Tkewness"+c,
                type=int,
                default=2**30)
        parser.add_argument("--init-lim",
                help="Factor base bound"+c,
                default=2**26)
        parser.add_argument("--init-lpb",
                help="Large prime bound"+c,
                default=64)
        parser.add_argument("--init-mfb",
                help="Cofactor bound"+c,
                default=100)
        parser.add_argument("--init-ncurves",
                help="ECM effort in cofactorization"+c,
                default=80)
        parser.add_argument("--init-I",
                help="Sieving range"+c,
                default=14)
        parser.add_argument("--init-minB1",
                help="ECM first B1" + c,
                default=200)
        parser.add_argument("--init-mineff",
                help="ECM minimal effort" + c,
                default=1000)
        parser.add_argument("--init-maxeff",
                help="ECM maximal effort" + c,
                default=100000)
        parser.add_argument("--init-side",
                help="Side of the bootstrap (when there is no rational side)",
                default=1)
        # Slave las processes in the initial step.
        parser.add_argument("--slaves",
                help="Number of slaves to use",
                type=int, default=1)
        # In case we used an external process
        parser.add_argument("--external-init",
                help="Use precomputed external data for the descent bootstrap",
                type=str,
                default=None)

    def __init__(self, general, args):
        self.general = general
        self.logDB = general.logDB

        if args.external_init != None:
            self.external = args.external_init
            if not os.path.exists(self.external):
                raise NameError("Given external file for init does not exist")
        else:
            self.external = None
            self.tkewness = int(args.init_tkewness)
            self.lim      = int(args.init_lim)
            self.lpb      = int(args.init_lpb)
            self.mfb      = int(args.init_mfb)
            self.ncurves  = int(args.init_ncurves)
            self.I        = int(args.init_I)
            self.side     = int(args.init_side)
            self.mineff   = int(args.init_mineff)
            self.maxeff   = int(args.init_maxeff)
            self.minB1    = int(args.init_minB1)
            self.slaves   = int(args.slaves)
            # the final step needs to know the init side as well.
            general.init_side = int(args.init_side)

    def __isqrt(self, n):
        x = n
        y = (x + 1) // 2
        while y < x:
            x = y
            y = (x + n // x) // 2
        return x

    def __myxgcd(self, a, b, T):
        assert type(a) == int
        assert type(b) == int
        assert type(T) == int
        ainit = a
        binit = b
        bound = self.__isqrt(b*T)
        x = 0
        lastx = 1
        y = 1
        lasty = 0
        while abs(b) > bound:
            q = a // b
            r = a % b
            a = b
            b = r
            newx = lastx - q*x
            lastx = x
            x = newx
            newy = lasty - q*y
            lasty = y
            y = newy
        return [ [ b, x ], [ a, lastx ] ]

    def use_external_data(self, z):
        fil = open(self.external, "r")
        rrr = fil.read()
        fil.close()
        lines = rrr.splitlines()
        e = int(lines[0])
        Num = int(lines[1])
        Den = int(lines[2])
        ## check that we are talking about the same z!
        p = general.p()
        zz = pow(z, e, p)
        assert (zz*Den-Num) % p == 0
        general.initrandomizer = e       # for later use
        fnum = [ int(x) for x in lines[3].split() ]
        fden = [ int(x) for x in lines[4].split() ]
        large_q = [ int(x) for x in lines[5].split() ]
        descrelfile = lines[6]

        ## create todolist from fnum and fden, skipping primes of the
        ## large_q list
        prefix = general.prefix() + ".descent.%s.init." % general.short_target()
        todofilename = os.path.join(general.datadir(), prefix + "todo")
        with open(todofilename, "w") as f:
            for q in fnum + fden:
                if q in large_q:
                    continue
                if self.logDB.has(q,-1,0):
                    continue
                logq = math.ceil(math.log(q, 2))
                print("Will do further descent for %d-bit rational prime %d"
                        % (logq, q))
                # las can understand when the rational root is missing
                f.write("0 %d\n" % q)
        fil = open(descrelfile, "r")
        rrr = fil.read()
        fil.close()
        lines = rrr.splitlines()
        with open(todofilename, "a") as f:
            for line in lines:
                foo = re.match("^Taken: ([0-9\-]+),([0-9\-]+):([0-9a-fA-F,]+):([0-9a-fA-F,]+)", line)
                assert foo
                foog = foo.groups()
                a = int(foog[0])
                b = int(foog[1])
                list_p = [[int(x, 16) for x in foog[i].split(",") ] for i in [2, 3]]
                for side in range(2):
                    for p in list_p[side]:
                        if p in large_q:
                            continue
                        if side == 0:
                            if not self.logDB.has(p,-1,0):
                                f.write("0 %d\n" % p)
                        else:
                            if b % p == 0:
                                continue
                            ideal = ideals_above_p(p, 1, a, b, side, general)
                            if ideal.get_log() != None:
                                continue
                            else:
                                r = a_over_b_mod_p(a, b, p)
                                f.write("1 %d %d\n" % (p, r))
        return todofilename, [Num, Den, fnum, fden], descrelfile

    def do_descent_for_real(self, z, seed):
        p = general.p()
        bound = p.bit_length() // 2 + 20
        # make the randomness deterministic to be able to replay
        # interrupted computations.
        random.seed(seed)
        general.initrandomizer = random.randrange(p)
        while True:
            zz = pow(z, general.initrandomizer, p)
            gg = self.__myxgcd(zz, p, self.tkewness)
            if (gg[0][0].bit_length() < bound and
                    gg[1][0].bit_length() < bound and
                    gg[0][1].bit_length() < bound and
                    gg[1][1].bit_length() < bound):
                break
            print ("Skewed reconstruction. Let's randomize the input.")
            general.initrandomizer = random.randrange(p)

        tmpdir = general.tmpdir()
        prefix = general.prefix() + ".descent.%s.%s.init." % (general.short_target(), seed)

        polyfilename = os.path.join(tmpdir, prefix + "poly")
        with open(polyfilename, 'w') as f:
            f.write("n: %d\n" % p)
            f.write("skew: 1\n")
            f.write("c1: %d\n" % gg[0][0])
            f.write("c0: %d\n" % gg[1][0])
            f.write("Y1: %d\n" % gg[0][1])
            f.write("Y0: %d\n" % gg[1][1])

        print ("--- Sieving (initial) ---")

        def relation_filter(data):
            line = data.strip()
            if re.match(r"^[^#]-?\d+,\d+:(\w+(,\w+)*)?:(\w+(,\w+)*)?$", line):
                return line

        relsfilename = os.path.join(general.datadir(), prefix + "rels")

        if os.path.exists(relsfilename):
            sources = [ (relsfilename, []) ]
        else:
            fbcfilename = os.path.join(tmpdir, prefix + "fbc")
            call_common = [ general.las_bin(),
                      "-poly", polyfilename,
                      "-lim0", self.lim,
                      "-lim1", self.lim,
                      "-lpb0", self.lpb,
                      "-lpb1", self.lpb,
                      "-mfb0", self.mfb,
                      "-mfb1", self.mfb,
                      "-ncurves0", self.ncurves,
                      "-ncurves1", self.ncurves,
                      "-fbc", fbcfilename,
                      "-I", self.I
                ]
            def fbc_call():
                call_that = call_common + [
                        "-q0", self.tkewness,
                        "-q1", self.tkewness,
                        "-nq", 0,
                        "-t", "machine,1,pu" if has_hwloc else "4"
                ]
                call_that = [str(x) for x in call_that]
                return call_that

            def construct_call(q0,q1):
                call_that = call_common + [
                          "-q0", q0,
                          "-q1", q1,
                          "--exit-early", 2,
                          "-t", "auto" if has_hwloc else "4"
                ]
                call_that = [str(x) for x in call_that]
                return call_that

            call_params = [(os.path.join(relsfilename+"."+str(i)), # outfile
                            self.tkewness+100000*i, #q0
                            self.tkewness+100000*(i+1)) for i in range(self.slaves)] #q1

            if not os.path.exists(fbcfilename):
                all_ok = True
                for t in call_params:
                    if not os.path.exists(t[0]):
                        all_ok = False
                        break

                if all_ok:
                    print (" - Using %s" % fbcfilename)
                else:
                    print (" - Factor base cache -")
                    with open(os.devnull, 'w') as devnull:
                        subprocess.check_call(fbc_call(), stdout=devnull)
                    print (" - done -")

            # Whether or not the output files are already present, this
            # will do the right thing and run the new processes only if
            # needed.
            sources = [(outfile, construct_call(q0,q1)) for (outfile,q0,q1) in call_params]

        rel_holder = object_holder()

        def consume(rel_holder, idx, line):
            if line[0] != '#':
                rel_holder.set((idx, line))
                return True
            if re.match("^# (Now sieving.*q=|\d+ relation)", line):
                sys.stdout.write('\n')
                print(line.rstrip())
                sys.stdout.flush()

        monitor_important_files(sources,
                                consume,
                                (rel_holder,),
                                temporary_is_reusable=True)
                       
        sys.stdout.write('\n')
        if rel_holder.v is None:
            print("No relation found!")
            print("Trying again with another random seed...")
            return None, None, None

        idx,rel = rel_holder.v

        if not os.path.exists(relsfilename):
            shutil.copyfile(sources[idx][0],relsfilename)
        
        print("Taking relation %s\n" % rel)
        rel = rel.split(':')
        a,b = [int(x) for x in rel[0].split(',')]

        Num = a*gg[0][0] + b*gg[1][0]
        Den = a*gg[0][1] + b*gg[1][1]
        assert (zz*Den-Num) % p == 0

        factNum = [ int(x, 16) for x in rel[2].split(',') ]
        factDen = [ int(x, 16) for x in rel[1].split(',') ]
        print(Num, Den, factNum, factDen)

        assert(abs(Num) == functools.reduce(lambda x,y:x*y,factNum,1))
        assert(abs(Den) == functools.reduce(lambda x,y:x*y,factDen,1))

        lc_ratpol = int(general.poly_data()["Y"][1])
        for q in factNum + factDen:
            if not self.logDB.has(q,-1,0):
                if lc_ratpol % q == 0:
                    print("Would need to descend %s which divides the lc of the rational poly." % q)
                    print("Trying again with a new seed.")
                    return None, None, None

        todofilename = os.path.join(general.datadir(), prefix + "todo")

        if not os.path.exists(todofilename):
            with open(todofilename, "w") as f:
                for q in factNum + factDen:
                    if self.logDB.has(q,-1,0):
                        continue
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent for %d-bit rational prime %d"
                            % (logq, q))
                    # las can understand when the rational root is missing
                    f.write("0 %d\n" % q)
        else:
            with open(todofilename, "r") as f:
                for line in f:
                    side,q = line.strip().split(' ')
                    q=int(q)
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent for %d-bit rational prime %d" % (logq, q))


        return todofilename, [Num, Den, factNum, factDen], None

    def do_descent_nonlinear(self, z):
        p = general.p()
        tmpdir = general.tmpdir()
        prefix = general.prefix() + ".descent.%s.upper." % general.short_target()
        polyfilename = os.path.join(tmpdir, prefix + "poly")
        if general.extdeg() == 1:
            zz = [ z ]
        else:
            zz = z
        call_that = [ general.descentinit_bin(),
                "-poly", general.poly(),
                "-mt", 4,
                "-minB1", self.minB1,
                "-mineff", self.mineff,
                "-maxeff", self.maxeff,
                "-side", self.side,
                "-extdeg", general.extdeg(),
                "-lpb", self.lpb,
                "-seed", 42,
                "-jl",
                p ] + zz
        call_that = [str(x) for x in call_that]
        initfilename = os.path.join(general.datadir(), prefix + "init")

        has_winner = object_holder(True)

        def consume(has_winner, general, idx, line):
            foo = re.match("^Youpi: e = (\d+) is a winner", line)
            if foo:
                has_winner.set(True)
                general.initrandomizer = int(foo.groups()[0])
            foo = re.match("^U = ([0-9\-,]+)", line)
            if foo:
                general.initU = [ int(x) for x in foo.groups()[0].split(',') ]
            foo = re.match("^V = ([0-9\-,]+)", line)
            if foo:
                general.initV = [ int(x) for x in foo.groups()[0].split(',') ]
            foo = re.match("^u = ([0-9]+)", line)
            if foo:
                general.initu = int(foo.groups()[0])
            foo = re.match("^v = ([0-9]+)", line)
            if foo:
                general.initv = int(foo.groups()[0])
            foo = re.match("^fac_u = ([, 0-9]+)", line)
            if foo:
                general.initfacu = [ [ int(y) for y in x.split(',') ] for x in foo.groups()[0].split(' ') ]
            foo = re.match("^fac_v = ([, 0-9]+)", line)
            if foo:
                general.initfacv = [ [ int(y) for y in x.split(',') ] for x in foo.groups()[0].split(' ') ]

        monitor_important_files(
                [(initfilename, call_that)],
                consume,
                (has_winner, general))

        if not has_winner.v:
            raise ValueError("initial descent failed for target %s" % zz)

        todofilename = os.path.join(general.datadir(), prefix + "todo")
        print(general.initfacu)
        print(general.initfacv)

        if not os.path.exists(todofilename):
            with open(todofilename, "w") as f:
                for ideal in general.initfacu + general.initfacv:
                    q = ideal[0]
                    r = ideal[1]
                    if self.logDB.has(q,r,self.side):
                        continue
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent for %d-bit prime %d"
                            % (logq, q))
                    f.write("%d %d %d\n" % (self.side, q, r))
        else:
            with open(todofilename, "r") as f:
                for line in f:
                    ll = line.strip().split(' ')
                    side = ll[0]
                    q = ll[1]
                    q=int(q)
                    logq = math.ceil(math.log(q, 2))
                    print("Will do further descent for %d-bit prime %d" % (logq, q))


        return todofilename, [general.initU, general.initV,
                general.initfacu, general.initfacv], None

    def do_descent(self, z):
        if not self.external:
            if general.has_rational_side():
                seed=42
                while True:
                    tdf, spl, frf = self.do_descent_for_real(z, seed)
                    if tdf != None:
                        return tdf, spl, frf
                    else:
                        seed += 1
            else:
                return self.do_descent_nonlinear(z)
        else:
            return self.use_external_data(z)

class DescentMiddleClass(object):
    def declare_args(parser):
        # TODO: import default values from the sieving parameters.
        parser.add_argument("--descent-hint",
            help="Hintfile for the descent",
            required=True,  # TODO: fall back on default values.
            )
        parser.add_argument("--I",
            help="Default value for I (must match hint file)",
            required=True,
            type=int)
        for side in range(2):
            parser.add_argument("--mfb%d" % side,
                    help="Default cofactor bound on side %d" % side,
                    required=True,
                    type=int)
            parser.add_argument("--lim%d" % side,
                    help="Default factor base bound on side %d (must match hint file)" % side,
                    required=True,
                    type=int)

    def __init__(self, general, args):
        self.general = general
        self.args = args
        # We need to do some safety checking
        values_I=set()
        values_lim0=set()
        values_lim1=set()
        values_I.add(args.I)
        values_lim0.add(args.lim0)
        values_lim1.add(args.lim1)
        with open(args.descent_hint, 'r') as file:
            for line in file:
                if re.match("^\s*#", line):
                    continue
                if re.match("^\s*$", line):
                    continue
                line = line.strip()
                foo = re.match("^.*I=(\d+)\s+(\d+),[\d.,]+\s+(\d+),[\d.,]+$",
                        line)
                if not foo:
                    print("Warning, parse error in hint file at line:\n" + line)
                    continue
                I,lim0,lim1 = foo.groups()
                values_I.add(int(I))
                values_lim0.add(int(lim0))
                values_lim1.add(int(lim1))
        if len(values_lim0)>1:
            raise ValueError("lim0 values should match between cmdline and hint file")
        if len(values_lim1)>1:
            raise ValueError("lim1 values should match between cmdline and hint file")
        if len(values_I)>1:
            raise ValueError("I values should match between cmdline and hint file")
        print("Consistency check for las_descent passed")
        print("\tI=%d" % values_I.pop())
        print("\tlim0=%d" % values_lim0.pop())
        print("\tlim1=%d" % values_lim1.pop())


    def do_descent(self, todofile):
        tmpdir = general.tmpdir()
        prefix = general.prefix() + ".descent.%s.middle." % general.short_target()

        f = open(todofile, 'r')
        ntodo = len(list(f))
        f.close()
        print ("--- Sieving (middle, %d rational primes) ---" % ntodo)
        s=general.lasMiddle_base_args()
        if args.descent_hint:
            s += [ "--descent-hint-table", args.descent_hint ]
        s += [
                "--I", self.args.I,
                "--lim0", self.args.lim0,
                "--lim1", self.args.lim1,
                "--lpb0", general.lpb0(),
                "--mfb0", self.args.mfb0,
                "--lpb1", general.lpb1(),
                "--mfb1", self.args.mfb1,
                "-t", "machine,1,pu" if has_hwloc else "4"
             ]
        s += [ "--todo", todofile ]
        call_that=[str(x) for x in s]
        relsfilename = os.path.join(general.datadir(), prefix + "rels")

        printing = object_holder(False)
        failed = []
        
        def consume(printing, failed, idx, line):
            if re.match("^# taking path", line):
                print(line.rstrip())
            elif re.match("^# END TREE", line):
                print("")
                printing.unset()
            elif printing.v:
                print(line.rstrip())
                foo = re.match("# FAILED (\d+\@\d+)", line)
                if foo:
                    failed.append(foo.groups()[0])
            elif re.match("^# BEGIN TREE", line):
                print("")
                printing.set(True)

        monitor_important_files(
                [(relsfilename, call_that)],
                consume,
                (printing, failed)
                )

        if failed:
            raise RuntimeError("Failed descents for: " + ", ".join(failed))

        return relsfilename

def prime_ideal_mixedprint(pr):
    p = pr[0]
    side = pr[1]
    if side == 0:
        machine = "%x 0 rat" % p
        human = "0,%d" % p
    else:
        r = pr[2]
        machine = "%x %d %x" % pr
        human = "%d,%d,%d" % (side,p,r)
    return machine,human


class DescentLowerClass(object):
    def declare_args(parser):
        pass
    def __init__(self, general, args):
        self.general = general
        self.args = args

    def __count_multiplicites(self, L):
        LL = []
        prev_p = L[0]
        m = 1
        for i in range(1,len(L)):
            p = L[i]
            if p == prev_p:
                m += 1
            else:
                LL.append([prev_p, m])
                prev_p = p
                m = 1
        LL.append([prev_p, m])
        return LL

    def do_descent(self, relsfile, initial_split):
        args = parser.parse_args()
        tmpdir = general.tmpdir()
        prefix = general.prefix() + ".descent.%s.lower." % general.short_target()
        relsforSM = os.path.join(tmpdir, prefix + "relsforSM")
        SMfile = os.path.join(tmpdir, prefix + "SM")

        # Read descent relations
        descrels = []
        for rfile in relsfile:
            with open(rfile, 'r') as file:
                with open(relsforSM, 'a') as fileSM:
                    for line in file:
                        foo = re.match("^Taken: (-?\d+),(-?\d+):", line)
                        if foo:
                            r = line.split(':')[1:]
                            r[0] = r[0].lstrip()
                            fileSM.write(r[0] + ":" + r[1] + ":" + r[2])
                            a,b = r[0].split(',')
                            a=int(a)
                            b=int(b)
                            list_p = [ [], [] ]
                            for side in range(2):
                                for p in r[side+1].strip().split(','):
                                    list_p[side].append(int(p, 16))
                            list_p = [ self.__count_multiplicites(list_p[0]),
                                    self.__count_multiplicites(list_p[1])]
                            descrels.append(([a,b], list_p))
        nrels = len(descrels)
        print ("--- Final reconstruction (from %d relations) ---" % nrels)

        # Compute SM
        call_that = [ general.sm_simple_bin(),
                        "-poly", general.poly(),
                        "-inp", relsforSM,
                        "-out", SMfile,
                        "-ell", general.ell()
                    ]
        if self.args.sm_mode is not None:
            call_that += [ "-sm-mode", self.args.sm_mode ]
        call_that = [str(x) for x in call_that]
        print("command line:\n" + " ".join(call_that))
        with open(os.devnull, 'w') as devnull:
            subprocess.check_call(call_that, stderr=devnull)

        SM = []
        with open(SMfile, 'r') as file:
            for line in file:
                r = line.split()
                sm = [ int(x) for x in r ]
                SM.append(sm)
        assert len(SM) == nrels

        # Reverse the order of relations to get only one unknown log
        # per relation while processing them
        descrels.reverse()
        SM.reverse()

        # Fill-in the log database
        logDB = general.logDB
        irel = 0
        for rel in descrels:
            unk = None
            a, b = rel[0]
            list_p = rel[1]
            acc_log = 0
            if logDB.fullcolumn != None:
                acc_log += logDB.fullcolumn
            sm = SM[irel]
            for i in range(len(sm)):
                acc_log += logDB.allSM(i)*sm[i]
            for side in range(2):
                for p, k in list_p[side]:
                    ideal = ideals_above_p(p, k, a, b, side, general)
                    log = ideal.get_log()
                    if log == None:
                        if unk != None:
                            raise ValueError(
                        "Two unknown ideals in relation a,b=%d,%d: %d (side %d) and %d (side %d)"
                        % (a, b, unk[0], unk[2], p, side))
                        else:
                            unk = [p, a_over_b_mod_p(a, b, p), side]
                    else:
                        acc_log += ideal.get_log()
            acc_log = acc_log % general.ell()
            if unk == None:
                assert acc_log == 0
            else:
                log = general.ell() - acc_log
                print ("Deduced log of (%d, %d, %d) from rel: %d"
                        % (unk[0], unk[1], unk[2], log))
                logDB.add_log(unk[0], unk[1], unk[2], log)
            irel += 1

        if general.has_rational_side():
            ## Deduce the log of the target
            Num, Den, factNum, factDen = initial_split
            log_target = 0
            errors=[]
            for p in factNum:
                lp = logDB.get_log(p, -1, 0)
                if lp is None:
                    errors.append(p)
                else:
                    log_target = log_target + lp
            for p in factDen:
                lp = logDB.get_log(p, -1, 0)
                if lp is None:
                    errors.append(p)
                else:
                    log_target = log_target - lp
            if len(errors):
                msg = "Some logarithms missing:\n"
                msg += "\n".join(["\t"+str(x) for x in errors])
                raise RuntimeError(msg)
            p=general.p()
            ell=general.ell()
            log_target = log_target % ell
            if general.initrandomizer != 1:
                # divide result by randomizer modulo ell
                multiplier = pow(general.initrandomizer, ell-2, ell)
                log_target = (log_target * multiplier) % ell
            print("# p=%d" % p)
            print("# ell=%d" % ell)
            print("log(2)=%d" % logDB.get_log(2, -1, 0))
            print("log(3)=%d" % logDB.get_log(3, -1, 0))
            print("# target=%s" % args.target)
            print("log(target)=%d" % log_target)
            check_result(2, logDB.get_log(2, -1, 0), int(args.target), log_target, p, ell)
        else:
            ## No rational side; more complicated.
            # We need to compute the SMs for U and V.
            polyforSM = os.path.join(tmpdir, prefix + "polyforSM")
            SM2file = os.path.join(tmpdir, prefix + "SM2")

            with open(polyforSM, 'w') as f:
                for poly in [ general.initU, general.initV ]:
                    f.write("p %d" % (len(poly)-1))
                    for c in poly:
                        f.write(" %d" % c)
                    f.write("\n")

            call_that = [ general.sm_simple_bin(),
                            "-poly", general.poly(),
                            "-inp", polyforSM,
                            "-out", SM2file,
                            "-ell", general.ell()
                        ]
            call_that = [str(x) for x in call_that]
            print("command line:\n" + " ".join(call_that))
            with open(os.devnull, 'w') as devnull:
                subprocess.check_call(call_that, stderr=devnull)

            SM2 = []
            with open(SM2file, 'r') as file:
                for line in file:
                    r = line.split()
                    sm = [ int(x) for x in r ]
                    SM2.append(sm)
            assert len(SM2) == 2

            ell = general.ell()
            vlog = [0, 0]
            factored = [ general.initfacu, general.initfacv ]
            for i in range(0,2):
                for xx in factored[i]:
                    vlog[i] += logDB.get_log(xx[0], xx[1], general.init_side)
                ind_shift = 0
                if general.init_side == 1:
                    ind_shift = logDB.nSM(0)
                for j in range(logDB.nSM(general.init_side)):
                    vlog[i] += logDB.SM(general.init_side,j)*SM2[i][ind_shift+j]
                vlog[i] = vlog[i] % ell

            log_target = (vlog[0] - vlog[1]) % ell
            multiplier = pow(general.initrandomizer, ell-2, ell)
            log_target = (log_target * multiplier) % ell
            print("# p=%d" % general.p())
            print("# ell=%d" % ell)
            print("# target=%s" % args.target)
            print("log(target)=%d" % log_target)



# http://stackoverflow.com/questions/107705/disable-output-buffering
# shebang takes only one arg...
# python3 doesn't grok sys.stdout = os.fdopen(sys.stdout.fileno(), 'w', 0)
# setting PYTHONUNBUFFERED here is too late.
class drive_me_crazy(object):
    def __init__(self, stream, timestamp=False):
        self.stream = stream
        self.eol = 1
        self.timestamp = timestamp
    def write(self, data):
        if self.timestamp:
            p=0
            while len(data) > p:
                d = data.find('\n', p)
                if d < 0:
                    break
                if self.eol:
                    self.stream.write(time.asctime() + " ")
                self.stream.write(data[p:d+1])
                self.eol=True
                p = d + 1
            if len(data) > p:
                if self.eol:
                    self.stream.write(time.asctime() + " ")
                self.stream.write(data[p:])
                self.eol= False
        else:
            self.stream.write(data)
        self.stream.flush()
    def __getattr__(self, attr):
        return getattr(self.stream, attr)


if __name__ == '__main__':

    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Descent initialization for DLP")

    # Required
    parser.add_argument("--target", help="Element whose DL is wanted",
            type=str, required=True)
    parser.add_argument("--timestamp",
            help="Prefix all lines with a time stamp",
            action="store_true")


    GeneralClass.declare_args(parser)
    DescentUpperClass.declare_args(parser)
    DescentMiddleClass.declare_args(parser)
    DescentLowerClass.declare_args(parser)

    args = parser.parse_args()

    sys.stdout = drive_me_crazy(sys.stdout, args.timestamp)

    las_bin = os.path.join(args.cadobindir, "sieve", "las")
    cp = subprocess.Popen([ las_bin, "-help" ],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE)
    if re.search("unused, needs hwloc", cp.stderr.read().decode()):
        has_hwloc = False
    else:
        has_hwloc = True

    general = GeneralClass(args)

    if general.target() == 1:
        # the re-randomization does not work for target=1
        print("# p=%d" % general.p())
        print("# ell=%d" % general.ell())
        print("# target=%s" % args.target)
        print("log(target)=0")
    else:
        init = DescentUpperClass(general, args)
        middle = DescentMiddleClass(general, args)
        lower = DescentLowerClass(general, args)

        todofile, initial_split, firstrelsfile = init.do_descent(general.target())
        relsfile = middle.do_descent(todofile)
        if firstrelsfile:
            lower.do_descent([firstrelsfile, relsfile], initial_split)
        else:
            lower.do_descent([relsfile], initial_split)

    general.cleanup()
