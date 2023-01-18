# read a matrix output by replay (in txt mode)
# M=read_matrix("c60.sparse.txt")
# L=read_matrix("c60.sparseL.txt")
# R=read_matrix("c60.sparseR.txt")
def read_matrix(f):
   f = open(f,"r")
   s = f.readline()
   s = s.split()
   assert len(s) == 2, "len(s) == 2"
   nrows = ZZ(s[0])
   ncols = ZZ(s[1])
   M = [[nrows,ncols]]
   maxcol = -1
   J = []
   w = 0
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split()
      l = [ZZ(x) for x in s]
      k = l[0] # k should be the number of elements in tail(l)
      assert len(l) == k+1, "len(l) == k+1"
      l = l[1:]
      l.sort()
      M.append(l)
      w += len(l)
      if l!=[]:
        maxcol = max(maxcol, max(l))
      for j in l:
         if not j in J:
            J.append(j)
   assert maxcol<ncols, "maxcol<ncols"
   print ("dimensions:", M[0])
   print ("number of different column elements:", len(J))
   print ("weight:", w)
   f.close()
   return M

# MM=sparse_mul(L,R)
def sparse_mul(L,R):
   n = L[0][0]
   m = L[0][1]
   j = R[0][0]
   assert m==j, "m==j"
   p = R[0][1]
   M = dict()
   L = L[1:]
   R = R[1:]
   assert len(L)==n, "len(L)==n"
   assert len(R)==m, "len(R)==m"
   for i in range(n):
      for j in L[i]:
         for k in R[j]:
            # multiply L[i,j] by R[j,k] and add to M[i,k]
            if not (i,k) in M.keys():
               M[i,k] = 0
            M[i,k] = (M[i,k]+1) % 2
   # remove 0 values
   Mcopy = dict()
   w = 0
   for (i,k) in M.keys():
      if M[i,k]!=0:
         Mcopy[i,k] = M[i,k]
         w += 1
   print ("weight=", w)
   M = [[n,p]]
   for i in range(n):
      M.append([])
   for (i,k) in Mcopy.keys():
      M[i+1].append(k)
   for i in range(n):
      M[i+1].sort()
   return M

# cmp_matrix(M,MM)
def cmp_matrix(M1,M2):
   if M1[0] != M2[0]:
      print ("dimensions differ:", M1[0], M2[0])
      return
   n = M1[0][0]
   diff = 0
   for i in range(n):
      if M1[i+1] != M2[i+1]:
         if diff == 0:
          print ("rows ", i, "differ")
         diff += 1
   print ("number of different rows:", diff)

# read index file
# I=read_index("c60.index")
def read_index(f):
   f = open(f,"r")
   s = f.readline()
   n = ZZ(s)
   print ("nrows:", n)
   M = [[n]]
   w = 0
   while true:
      s = f.readline()
      if s=='':
         break
      s = s.split()
      k = ZZ(s[0])
      s = s[1:]
      assert len(s) == k, "len(s) == k"
      l = []
      for c in s:
         l.append(ZZ('0x'+c))
      M.append(l)
      w += len(l)
   f.close()
   print ("weight:", w)
   return M

# M = read_purged("c60.purged",skip=32)
def read_purged(f,skip=0):
   f = open(f,"r")
   s = f.readline().split()
   nrows = ZZ(s[1])
   ncols = ZZ(s[2])
   M = []
   while true:
      s = f.readline()
      if s=='':
         break
      s = s[:-1]
      s = s.split(':')[1]
      s = s.split(',')
      l = [ZZ('0x'+x) for x in s]
      # remove skipped ideals
      l = [x for x in l if x >= skip]
      l.sort()
      M.append(l)
   f.close()
   return M

# M = read_purged("c60.purged",skip=32)
# Mprime = merge(M,"c60.history")
def merge(M,history):
   M = [copy(l) for l in M]
   f = open(history,"r")
   N = len(M)
   while true:
      s = f.readline()
      if s=='':
         break
      if s[0]=='#':
         continue
      s = s[:-1]
      s = s.split(' ')
      l = []
      remove = true
      for x in s:
         if x[0]!='#':
            i = ZZ(x)
            if i<0:
               remove = false
            l.append(i)
      # A line i1 i2 ... ik means that row i1 is added to i2, ..., ik, and
      # row i1 is removed afterwards (where row 0 is the first line in
      # *.purged.gz).
      i1 = l[0]
      l1 = M[i1]
      for i2 in l[1:]:
         # i2 += i1
         M[i2] = add_row (l1, M[i2])
      if remove:
         M[i1] = []
         N -= 1
         if N in [18723,17509,17382,17376,17375,15167,13665,12611,11817,11203]:
            print (N, weight(M))
   f.close()
   return M

def add_row(l1,l2):
   l = []
   n1 = len(l1)
   n2 = len(l2)
   i1 = i2 = 0
   while i1 < n1 and i2 < n2:
      if l1[i1]<l2[i2]:
         l.append(l1[i1])
         i1 += 1
      elif l2[i2]<l1[i1]:
         l.append(l2[i2])
         i2 += 1
      else: # cancellation
         i1 += 1
         i2 += 1
   while i1 < n1:
      l.append(l1[i1])
      i1 += 1
   while i2 < n2:
      l.append(l2[i2])
      i2 += 1
   return l
   
def weight(M):
   return add(len(l) for l in M)

def nrows(M):
   n = 0
   for l in M:
      if l!=[]:
         n += 1
   return n
