# Notes on using klone

#### Overview: klone is one of the UW supercomputers in the hyak system.

We own 400 cores (10 nodes with 40 cores each).  We are allocated 1 TB of storage for each node, so 10 TB total.

ssh to `klone.hyak.uw.edu` to get to it

`/gscratch/macc` is our working directory, and I have created my own directory inside that: parker.

There appears to be a klone1 so I will use that to be explicit (klo alias on my mac) although this is probably not required.  There is no klone2.  In contrast, on mox there is a mox1 and mox2.

---

#### From David Darr: klone requires only modest changes to Linux-ifort_mox.mk, renamed Linux-ifort_klone.mk (in LiveOcean_roms/LO_ROMS/Compilers).  The only difference actually is that the two instances of:
```
NC_CONFIG ?= nc-config
```
are now:
```
NC_CONFIG ?= nf-config
```

---

#### Settings for .bashrc (I added these to my .bashrc on klone)
```
module load intel/oneAPI
LODIR=/gscratch/macc/local
#OMPI=${LODIR}/openmpi-ifort
NFDIR=${LODIR}/netcdf-ifort
NCDIR=${LODIR}//netcdf-icc
export LD_LIBRARY_PATH=${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}
export PATH=/gscratch/macc/local/netcdf-ifort/bin:$PATH
export PATH=/gscratch/macc/local/netcdf-icc/bin:$PATH
#export PATH=/gscratch/macc/local/openmpi-ifort/bin:$PATH
```

---

#### Steps to compile:
```
cd LiveOcean_roms/LO_ROMS
srun -p compute -A macc --pty bash -l
make clean
make -f /gscratch/macc/parker/LiveOcean_roms/makefiles/[ex_name]/makefile
```
Then `logout` to get back to the usual shell.  You have to do this because the `srun` command logged you onto one of the compute nodes.

Note: `compute` in the srun command was `macc` on mox

---

#### Set up ssh-keygen to apogee

Log onto klone1 and do:
```
ssh-keygen
```
and hit return for most everything.  However, you may encounter a prompt like this:
```
Enter file in which to save the key (/mmfs1/home/pmacc/.ssh/id_rsa):
/mmfs1/home/pmacc/.ssh/id_rsa already exists.
Overwrite (y/n)?
```
Looking [HERE](https://www.hostdime.com/kb/hd/linux-server/the-guide-to-generating-and-uploading-ssh-keys), I found out that id_rsa is the default name that it looks for automatically.  You can name the key anything and then just refer to it when using ssh and etc. like:
```
ssh parker@apogee.ocean.washington.edu -i /path/to/ssh/key
```

So, since I want apogee to be the main user of klone I will chose to **overwrite** in the above.  This means that the old connection to boiler will no longer work.  When I did this it asked for a passphrase and I hit return.

Then I did:
```
ssh-copy-id parker@apogee.ocean.washington.edu
```
(it asks for my password)

And now I can ssh from klone to apogee without a password, and on apogee it had added a key with pmacc@klone1.hyak.local at the end to my ~/.ssh/authorized_keys.

On klone there is now an entry in ~/.ssh/known_hosts for apogee.ocean.washington.edu.

So, in summary: for going from klone1 to apogee it added to:
- ~/.ssh/known_hosts on klone (boiler and mox1 are also there), and
- ~/.ssh/authorized_keys on apogee

Question: if I now want to also allow other things like klone->perigee or mox->apogee, can I avoid overwriting my id_rsa?  Currently I rely on mox->boiler for the daily forecast.
