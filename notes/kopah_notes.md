# Notes on using the kopah file storage system

### kopah is part of the UW hyak computing system that also supports klone. It provides secure and accessible file storage using S3-API compatible object storage.

We (Parker MacCready and the Coastal Modeling Group) purchased a "Large Condo" of this resource in early 2025. This give us 300 TB of storage. The up-front cost was $14,141. There is also a monthly fee ($2.60/TB * 300 TB = $780/month) which amounts to $9,360/year. The contract is up to 6 years. Presumably after that you would have to pay a new up-front cost. Currently (1/2026) this is paid automatically on our NANOOS worktag. None of these expenses have IDC.

More kopah info: https://hyak.uw.edu/docs/storage/kopah

The goal of our using this resource is to transition away from maintaining our own servers (apogee and perigee) and instead do all our computing and storage on klone+kopah (or our laptops). This has three advantages. (1) We don't have to maintain the hardware. (2) We don't have to pay a system administrator such as David Darr, who has retired. (3) It is easy to make any file on kopah publically accessible using a URL, and this is much better for intereacting with stakeholders. The cost ($9k/year) is small compared to what we paid for a system administrator (~$35k/year = 2 months).

This resource is already being used in the LiveOcean system to make forecast extractions available to stakeholders (like NANOOS NVS) and for public access to monthly means: https://faculty.washington.edu/pmacc/LO/data_access.html. I also use it for making large custom extractions available, as a replacement for Dropbox or Google Drive. I use the tool LO/misc/kopah_share5.sh to do this, but note that this is hard-wired to use the "pm-share" bucket, so you would need to make your own version.

### A bit more about s3 buckets

Here is an AI summary of the bucket concept:

"Amazon S3 uses buckets as fundamental, top-level containers to organize data (objects) in the cloud, providing logical separation, security boundaries, and a way to apply consistent settings like access policies, versioning, and lifecycle rules, acting as the primary namespace for your data in AWS. Buckets are essential for managing vast amounts of unstructured data, offering scalability and unique global naming, with "folders" inside them actually being key prefixes for organizing objects."

You need to limit the number of buckets you create, but inside a bucket you can have a large number of folders, subfolders, and files.

For our "macc" group I would like to suggest that users initially create a bucket called "liveocean-[your UW NetID]", so for example mine is liveocean-pmacc. Inside of this you can put all the folders you need such as LO_output and LO_roms.

Moving forward all the LO tools ASSUME that this naming structure exists. In LO/driver I have tried to indicate the tools that conform to this new standard with the prefix "K_", e.g. K_driver_roms00.py.

Working with files in s3 buckets takes some getting used to. You cannot just access them with a simple path like you would for files on apogee and perigee. Instead, the command line tools described below need to be woven into your code.

## Command Line Tools

Excellent info about these is at https://hyak.uw.edu/docs/storage/cli.

**NOTE: I believe the command line tools I describe below are independent of the Globus file management system. Nonetheless, you should definitely read through and incorporate the tools described here AND in LO/Notes/globus_notes.md.**

There are two main command line tools: s3cmd and s5cmd. They are pretty similar but **s5cmd can be significantly faster**, so whenever possible that is what I use. These are linux tools, similar to the nco operators like ncks.

Both s3cmd and s5cmd are installed on klone, and they are also installed in the loenv conda environment. There is some inconsistency in access. In my login shell on klone both are available, but when I use srun to move to a compute node inly s3cmd is available. In either case they are both available when using loenv is activated. You can use e.g. "conda list s5cmd" to see the actual installed version number which occasionally is displayed as v0.0.0-dev during the conda update operation but which is really something like v2.3.0.

---

### s5cmd

This pretty much does everything that s3cmd (described below) does except: There is no direct "info" command in s5cmd that provides detailed metadata (like ACLs or headers) similar to s3cmd info. s5cmd is primarily optimized for high-performance bulk operations and lacks a specific metadata inspection sub-command.

[describe commands]

---

### s3cmd

To use this system you need to add a file ~/.s3cfg to any machine you use (klone, your laptop, etc.). The file has these lines:
```
[default]
host_base = s3.kopah.orci.washington.edu
host_bucket = s3.kopah.orci.washington.edu/%(bucket)
use_https = True
# Login credentials
access_key = [get it from Parker or Kate]
secret_key = [get it from Parker or Kate]
```

s3cmd is already one of the things installed by loenv.yml, and it is already on klone (so you don't need to be in loenv to use it on klone). We needed version 2.4.0 or better.

#### Examples of usage (all from the linux command line)

Get summary info about how much is in all the current buckets used by the macc group:
```
s3cmd du -H
```
List all objects in all buckets with human-readable sizes
```
s3cmd la -rH
```
Get detailed info for an object (file), including URL.
```
s3cmd info s3://[bucket]/[folder]/[object]
```
**Note that this reports http://s3.kopah.orci.washington.edu/ in the URL when it should be https://s3.kopah.uw.edu/! (https NOT http)**

Create a bucket. I think you need to create the bucket before you start putting things in it, but I am not sure.
```
s3cmd mb s3://[bucket name]
```
Bucket naming rules:
- Bucket names must be between 3 (min) and 63 (max) characters long.
- Bucket names can consist only of lowercase letters, numbers, periods (.), and hyphens (-).
- Bucket names must begin and end with a letter or number.
- You should start by making a bucket called liveocean-[your UW NetID].

Copy a file to a bucket:
```
s3cmd put [filename] s3://[bucket name]/[folder name if you like]/
```
Copy the file and make it publically accessible:
```
s3cmd put --acl-public [filename] s3://[bucket name]/[folder]/
```
**The URL to get the file is of the form:**

https://s3.kopah.uw.edu/[bucket name]/[folder]/[filename]

I am not sure about the details of granting more limited access, e.g. with a password. I also am not clear if the access can be object-by-object or if it must apply to a whole bucket.

Make a bucket and all its contents publically accessible if you forgot to do it when you first made it:
```
s3cmd setacl s3://[bucket name] --acl-public --recursive
```
Remove a bucket and its contents (if you remove one of my buckets I will not be happy):
```
s3cmd -rb s3://[bucket name] --force --recursive
```
When I was getting started with kopah I created too many buckets (150+?) and it wouldn't allow me to make any more. You can't use wild cards to specify a bunch of buckets in s3cmd so I wrote LPM/kopah/kopah_bucket_remover.py to get rid of them with a python program.

---
Copy the contents of a directory to a folder:
```
s3cmd sync [folder (no / needed)] s3://[bucket name]/
```
This would create s3://[bucket name]/[folder]/[all files]

---

### NOTES

There is still a lot to learn about using this system and I encourage you to do your own exploration. It will be a lot of work to rework the LO tools for pre and post processing to work with this system instead of the familiar apogee/perigee invironments. There are also likely new ways to speed up file access. I wrote LO/misc/speed_test.py to experiment with different ways of making a daily average of one model field from hourly files.

Icechunk may be a promising tool, created by Ryan Abernathy who is a pioneer in combining cloud computing with earth science:
https://earthmover.io/blog/icechunk.






