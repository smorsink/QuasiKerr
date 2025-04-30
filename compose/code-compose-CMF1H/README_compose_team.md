# **CompOSE**


### INSTALL :

in the normal case type "make" or "make compose"

If you want to build the version running on the WEB APP, type "make compose_for_web"


## DEPLOY ON THE WEB
after the `git push orgin`
ssh into the web server :
```
cd /data2/Supercompose/ComposeDjango/Eos_Tables/code/
git pull origin
make clean
```

Now you have two options:

1.
	Add to the begin of th Makefile :

```makefile
USE_BOUNDCHECK = 0
USE_HDF5 = 0
```
and then

```bash
make compose_for_web
```

2.
In a terminal :
```bash
export USE_BOUNDCHECK = 0
export USE_HDF5 = 0
make compose_for_web
```

On **  composecalc ** the code will not updated up to the chron execution.
