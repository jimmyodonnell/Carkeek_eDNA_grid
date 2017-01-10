load_package <-function(package_name){
  if(!require(package_name, character.only = TRUE)){
    install.packages(package_name)
    #call require again because install.packages simply installs the package, it doesn't load it
    require(package_name, character.only = TRUE)
  }
}
