\name{install_pro}
\alias{install_pro}
\title{Install ChannelAttributionPro binary tailored to your OS/arch/R}

\description{
Interactively installs the \pkg{ChannelAttributionPro} package by contacting the
ChannelAttribution build service, which returns a platform-appropriate binary (or
source) package URL for your current operating system, CPU architecture, and R
version. At the prompt you can enter:
\itemize{
  \item a \emph{token} (recommended), or
  \item a work/university \emph{email} to request a token (the function sends the request and exits).
}
}

\usage{
install_pro()
}

\details{
When called, the function asks you to enter your ChannelAttributionPro \emph{token}.
If you don't have one, you can enter your work/university email; the function will
request a token for that address and return invisibly, so you can re-run it once the
token arrives.

Processing steps:
\enumerate{
  \item Detect your platform (\code{os}, \code{os_vers}, \code{arch}) and R minor version.
  \item Send a form-encoded request to the ChannelAttribution builder endpoint to obtain a package URL.
  \item Install the returned package via \code{\link[utils]{install.packages}} with \code{repos = NULL}.
  \item On exit, send a small status message to a monitoring endpoint with minimal text (no local files).
}

\strong{Dependencies:}
This function requires \pkg{curl} (for HTTPS) and \pkg{jsonlite} (for JSON parsing). If either is missing,
the function aborts with a friendly message.

\strong{Networking and privacy:}
The function performs outbound HTTPS requests to \url{https://app.channelattribution.io}, sending only the minimal
parameters needed to resolve the correct binary and (on exit) a small status string for diagnostics. No local files
or datasets are uploaded.

\strong{Interactive prompt:}
Input is collected via base R \code{readline()} (no masking). In some environments the input may be visible.
}

\value{
Invisibly returns \code{NULL}. Progress and outcome messages are printed.
}

\section{Errors and messages}{
\itemize{
  \item \emph{Invalid or expired token (HTTP 401)}: You will see a message inviting you to contact
  \email{info@channelattribution.io} or request a new token.
  \item \emph{Network/SSL issues}: Reported as a concise \code{network_or_ssl_error: ...}.
  \item \emph{Unexpected builder response}: The function prints diagnostic information (platform and R version)
  that you can forward to support.
}
}

\seealso{
\code{\link[utils]{install.packages}}
}

\examples{
\dontrun{
## Typical interactive use:
install_pro()

## If you don't have a token at prompt:
##   - enter your work/university email to request one
##   - check your inbox (including spam)
##   - re-run install_pro() with the token
}
}

\author{ChannelAttribution Team}
\encoding{UTF-8}
\keyword{utilities}


