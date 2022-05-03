# README for post

#### The post code is very similar to the extract code, but it is designed specifically to work on one forecast, making all the expected post-processing products that I send to other people and the movies for the LiveOcean website.

In general the output files will go in:
```
out_dir = Ldir['LOo'] / 'post' / Ldir['gtagex'] / ('f' + Ldir['date_string']) / Ldir['job']
```
and this is created by `driver_post1.py` or `post_argfun.py` (and this is one reason we pass the -job argument to each `post_main.py`).
