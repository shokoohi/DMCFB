.onAttach <- function(lib, pkg) {
    packageStartupMessage(
    "DMCFB package, Version ", utils::packageVersion(
    "DMCFB"
    ), ", Released ", utils::packageDescription(
    "DMCFB"
    )$Date, "\n", utils::packageDescription(
    "DMCFB"
    )$Description,
    "\nBugReports: ", utils::packageDescription("DMCFB")$BugReports
    )
}
