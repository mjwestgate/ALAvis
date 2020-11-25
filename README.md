# ALAvis
A shiny app for linked spatial and taxonomic navigation of data from the Atlas of Living Australia.

To use:
```
# install.packages("shiny") # install shiny (if not done already)
library(shiny)
runGitHub("mjwestgate/ALAvis")
```

<b>Note: This app is still in early development. Many aspects either do not function or have suboptimal navigability.</b> A key issue for users is that <code>ALA4R</code> requires users to submit a registered email address; but this app does not implement that code, meaning that ALA searches will not work. Therefore users should treat this as a proof of concept only.