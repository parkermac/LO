# Globus Notes

## Globus is a high-perfromance file transfer system that is now fully supported in the UW hyak system.

The big deal is that it comes with a web-based file management app that allows drag-and-drop transfers freely across:

- your laptop
- apogee and perigee
- klone
- kopah

It easily and reliably handles huge transfers (many terrabytes!). You simply initiate the transfer in the website, then leave, and it will send you an email when it is done.

---

## Setup

Parker, end eventually Kate, needs to request mapping of your UW NetID to the macc kopah storage account. We do this by sending an email to help@uw.edu with "Kopah + Globus" in the subject line.

Once this is done you have to do a bit more installation, but the documentation in the hyak system is excellent.

#### perigee or apogee (Notes from David Darr)

Here is how to install this. If you have questions consult https://docs.globus.org/globus-connect-personal/install/linux/ or ask me. I find the online documentation somewhat misleading in a couple of places. Otherwise this is pretty straight forward.

download and extract (I did this in my /dat1/parker directory on apogee, since we try not to put too much in ~):

$  wget https://downloads.globus.org/globus-connect-personal/linux/stable/globusconnectpersonal-latest.tgz
```
tar xzf globusconnectpersonal-latest.tgz
```
this will produce a versioned globusconnectpersonal directory

replace `x.y.z` in the line below with the version number you see
```
cd globusconnectpersonal-x.y.z
```
run this command and follow instructions to setup:
```
./globusconnectpersonal -setup --no-gui
```
start (in background): 
```
./globusconnectpersonal -start &
```
check status:
```
./globusconnectpersonal -status
```
The next step is to edit the file: /home/$USER/.globusonline/lta/config-paths (by default globus only sees /home/$USER).  Note that you might have to create this file (it doesn't always seem to be there by default).

Here, is an example of how this looks on perigee after I added my perigee data directories:

```
more config-paths
```
```
~/,0,1
/data1/darr/,0,1
/data2/darr/,0,1
```
Of course you will need to add your own info here, not David's!

Once you get this going you will need to get "Globus mapped to Kopah". It basically entails writing UW-IT and requesting this (the "Setup" step above).

See https://hyak.uw.edu/docs/storage/gui for details.


