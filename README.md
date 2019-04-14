### gexp

This is the development code of the R package **gexp**.
You should use it if you want to contribute to its development:
testing unreleased versions, fixing bugs, writing code, etc.

To download, check and build it do the following in a terminal emulator:

> git clone  git://github.com/ivanalaman/gexp.git

> or

> git clone https://ivanalaman@github.com/gexp.git

After to clone it, to check, build and install do the following:
> R CMD check gexp

> R CMD build gexp

> R CMD INSTALL gexp_X.X-X.tar.gz

Or, you can install using devtools package as:

> library(devtools)

> install_github('ivanalaman/gexp')
