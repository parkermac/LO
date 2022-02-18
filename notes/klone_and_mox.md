# Notes on using klone and mox

#### Overview: klone and mox are UW supercomputers in the hyak system.

Here are examples of aliases I have on my mac ~/.bash_profile (equivalent to ~/.bashrc on the linux machines) to quickly get to my machines
```
alias klo='ssh pmacc@klone1.hyak.uw.edu'
alias mox1='ssh pmacc@mox1.hyak.uw.edu'
alias mox2='ssh pmacc@mox2.hyak.uw.edu'
alias pgee='ssh parker@perigee.ocean.washington.edu'
alias agee='ssh parker@apogee.ocean.washington.edu'
```
Note: klone1 is the same as klone.  If you just ssh to mox you end up randomly at either mox1 or mox2, which are the same machine except that they keep separate crontabs (see below). I always ssh to mox1 to avoid confusion.

`/gscratch/macc` is our working directory on both mox and klone, and I have created my own directory inside that: parker, where the whole LO system lives.

When you have a job running on klone you can check on it using:
```
squeue -A macc
```
or on mox:
```
squeue -p macc
```
---

#### Getting resource info

`hyakstorage` will give info about storage on klone.  Use `hyakstorage --help` to get more info on command options. Not yet working on mox.

To check on our disk allocation on mox you can also look in the file `/gscratch/macc/usage_report.txt` although this will be phased out soon.

`hyakalloc` will give info on the nodes we own.

**mox**: we own 196 cores (7 nodes with 28 cores each), plus another 64 (2 nodes with 32 cores each).

**klone**: we own 400 cores (10 nodes with 40 cores each). We are allocated 1 TB of storage for each node, so 10 TB total.

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
I assume the // in NCDIR path is a typo.

---

#### Steps to compile:

On klone:

```
cd LiveOcean_roms/LO_ROMS
srun -p compute -A macc --pty bash -l
make clean
make -f /gscratch/macc/parker/LiveOcean_roms/makefiles/[ex_name]/makefile
```
Then `logout` to get back to the usual shell.  You have to do this because the `srun` command logged you onto one of the compute nodes.

On mox the steps are only slightly different. The `compute` in the srun command is `macc` on **mox**:
```
cd LiveOcean_roms/LO_ROMS
srun -p macc -A macc --pty bash -l
make clean
make -f /gscratch/macc/parker/LiveOcean_roms/makefiles/[ex_name]/makefile
```

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

In the interests of tidying up I will chose to **overwrite** in the above.  This means that the old connection to boiler will no longer work.  When I did this it asked for a passphrase and I hit return.

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

Now I can run `ssh-copy-id` again for other computers, without having to do the `ssh-keygen` step.

---

#### Running things by cron

See LO/driver/crontabs for curent versions.  These are discussed more in LO/README.md.
