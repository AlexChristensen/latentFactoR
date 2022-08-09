.onload <- function(libname, pkgname)
{library.dynam("latentFactoR",package=pkgname,lib.loc=libname)}

.onAttach <- function(libname, pkgname)
{
    msg <- styletext(styletext(paste("\nlatentFactoR (version ", packageVersion("latentFactoR"), ")\n", sep = ""), defaults = "underline"), defaults = "bold")
    msg <- paste(msg,'\nSee `?simulate_factors` to get started\n')
    msg <- paste(msg,"\nFor bugs and errors, submit an issue to <https://github.com/AlexChristensen/latentFactoR/issues>")

    packageStartupMessage(msg)
}
