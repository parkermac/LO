# Notes on using the kopah file storage system

### kopah is part of the UW hyak computing system that also supports klone. It provides secure and accessible file storage using S3-API compatible object storage.

We (Parker MacCready and the Coastal Modeling Group) purchased a "Large Condo" of this resource in early 2025. This give us 300 TB of storage. The up-front cost was $14,141. There is also a monthly fee ($2.60/TB * 300 TB = $780/month) which amounts to $9,360/year. The contract is up to 6 years. Presumably after that you would have to pay a new up-front cost. Currently (1/2026) this is paid automatically on our NANOOS worktag. None of these expenses have IDC.

More kopah info: https://hyak.uw.edu/docs/storage/kopah

The goal of our using this resource is to transition away from maintaining our own servers (apogee and perigee) and instead do all our computing and storage on klone+kopah (or our laptops). This has three advantages. (1) We don't have to maintain the hardware. (2) We don't have to pay a system administrator such as David Darr, who has retired. (3) It is easy to make any file on kopah publically accessible using a URL, and this is much better for intereacting with stakeholders. The cost ($9k/year) is small compared to what we paid for a system administrator (~$35k/year = 2 months).

This resource is already being used in the LiveOcean system to make forecast extractions available to stakeholders (like NANOOS NVS) and for public access to monthly means: https://faculty.washington.edu/pmacc/LO/data_access.html. I also use it for making large custom extractions available, as a replacement for Dropbox or Google Drive. I use the tool LO/misc/kopah_share.sh to do this, but note that this is hard-wired to use the "pm-share" bucket, so you would probably want to make your own version.

### How to use it

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
There are many ways to interact with kopah files. The one I have settled on is the s3cmd module (a linux command line tool). This is already one of the things installed by loenv.yml, and it is already on klone (so you don't need to be in loenv to use it on klone). We needed version 2.4.0 or better.

Then to do all the file things you usually do in linux you use an s3cmd command. One trick is that every file has to be in a "bucket" (like a directory) and for reasons I don't understand you need to limit the number of buckets you create. But inside a bucket you can have a large number of folders, subfolders, and files. Here is an AI summary of the bucket concept:

"Amazon S3 uses buckets as fundamental, top-level containers to organize data (objects) in the cloud, providing logical separation, security boundaries, and a way to apply consistent settings like access policies, versioning, and lifecycle rules, acting as the primary namespace for your data in AWS. Buckets are essential for managing vast amounts of unstructured data, offering scalability and unique global naming, with "folders" inside them actually being key prefixes for organizing objects."

### Examples of usage (all from the linux command line)

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
**Note that this reports http in the URL when it should be https!**

---
Create a bucket. I think you need to create the bucket before you start putting things in it, but I am not sure.
```
s3cmd mb s3://[bucket name]
```
Bucket naming rules:
- Bucket names must be between 3 (min) and 63 (max) characters long.
- Bucket names can consist only of lowercase letters, numbers, periods (.), and hyphens (-).
- Bucket names must begin and end with a letter or number.

---
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

---
Make a bucket and all its contents publically accessible if you forgot to do it when you first made it:
```
s3cmd setacl s3://[bucket name] --acl-public --recursive
```
---
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






