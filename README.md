### gerexp

This is the development code of the R package **gerexp**.
You should use it if you want to contribute to its development:
testing unreleased versions, fixing bugs, writing code, etc.

To download, check and build it do the following in a terminal emulator:

> git clone  git://github.com/ivanalaman/gerexp.git

> or

> git clone https://ivanalaman@github.com/gerexp.git

After to clone it, to check, build and install do the following:
> R CMD check gerexp

> R CMD build gerexp

> R CMD INSTALL gerexp_X.X-X.tar.gz

Or, you can install using devtools package as:

> library(devtools)

> install_github('ivanalaman/gerexp')


