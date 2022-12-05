# read a matrix output by replay (in txt mode)
# M=read_matrix("/tmp/c60.sparse.txt")
# L=read_matrix("/tmp/c60.sparseL.txt")
def read_matrix(f):
   f = open(f,"r")
   s = f.readline()
   s = s.split()
   assert len(s) == 2, "len(s) == 2"
   nrows = ZZ(s[0])
   ncols = ZZ(s[1])
   M = [[nrows,ncols]]
   maxcol = -1
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split()
      l = [ZZ(x) for x in s]
      M.append(l)
      maxcol = max(maxcol, max(l))
   print ("nrows=", nrows, "ncols=", ncols, "maxcol=", maxcol)
   return M

