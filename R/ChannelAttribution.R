# ChannelAttribution: Markov model for online multi-channel attribution
# Copyright (C) 2015 -   Davide Altomare and David Loris <https://channelattribution.io>

# ChannelAttribution is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ChannelAttribution is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with ChannelAttribution.  If not, see <http://www.gnu.org/licenses/>.

.v=packageVersion("ChannelAttribution")
 
.onAttach = function(libname, pkgname) {

 packageStartupMessage(paste0("ChannelAttribution ",.v))
 packageStartupMessage("*** Looking to run more advanced attribution? Install ChannelAttribution Pro for free running install_pro(). Visit https://channelattribution.io for more info.")

}

.message_pro="*** Install ChannelAttribution Pro for free running install_pro(). Visit https://channelattribution.io for more info. Set flg_pro=FALSE to hide this message."

heuristic_models=function(Data, var_path, var_conv, var_value=NULL, sep=">", flg_pro=TRUE){

 if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
  print("Data must be a data.frame or a data.table")
 } 
 
 if(is.character(var_path)){
	if(!var_path%in%names(Data)){
	 print("var_path must be a column of Data")
	}
 }else{
  print("var_path must be a string")
 }
 if(is.character(var_conv)){
	if(!var_conv%in%names(Data)){
	 print("var_conv must be a column of Data")
	}
 }else{
  print("var_conv must be a string")
 }
 
 if(!is.null(var_value)){
  if(!var_value%in%names(Data)){
   print("var_value must be a column of Data")
  }
 }

 if(length(sep)>1){stop("Separator must have length 1")}

 if(is.null(var_value)){var_value="0"}

 res=.Call("heuristic_models_cpp", Data, var_path, var_conv, var_value, sep)
 
 if(flg_pro==TRUE){
  print("*** Install ChannelAttribution Pro for free! Run install_pro(). Set flg_pro=FALSE to hide this message.")
 }
 
 return(as.data.frame(res)) 

}	

choose_order=function(Data, var_path, var_conv, var_null, max_order=10, sep=">", ncore=1, roc_npt=100, plot=TRUE, flg_pro=TRUE){
 
 if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
  print("Data must be a data.frame or a data.table")
 } 
 
 if(is.character(var_path)){
	if(!var_path%in%names(Data)){
	 print("var_path must be a column of Data")
	}
 }else{
  print("var_path must be a string")
 }
 
 if(is.character(var_conv)){
	if(!var_conv%in%names(Data)){
	 print("var_conv must be a column of Data")
	}
 }else{
  print("var_conv must be a string")
 }
  
 if(!is.null(var_null)){
  if(!var_null%in%names(Data)){
   print("var_null must be a column of Data")
  }
 }

 if(length(sep)>1){stop("sep must have length 1")}
 if(ncore<1){stop("ncore must be >= 1")}
 if(roc_npt<10){stop("roc_npt must be >= 10")}
 if(!plot%in%c(0,1)){stop("plot must be FALSE or TRUE")}

 res=.Call("choose_order_cpp", Data, var_path, var_conv, var_null, max_order, sep, ncore, roc_npt)
 
 ck=res$auc$order[res$auc$order!=0]
 res$auc$order=res$auc$order[ck]
 res$auc$auc=res$auc$auc[ck]
 res$auc$pauc=res$auc$pauc[ck]
 
 best_order=res$auc$order[res$auc$pauc==max(res$auc$pauc)]
 
 if(best_order==max_order){
  print(paste0("Suggested order not found. Try increasing max_order."))
 }else{
  print(paste0("Suggested order: ", res$auc$order[res$auc$pauc==max(res$auc$pauc)]))
 }
 
 if(plot=="TRUE"){
  plot(res$auc$order,res$auc$pauc,type="l",xlab="order",ylab="penalized auc",main="PENALIZED AUC")
 }
 
 res[['suggested_order']]=best_order
 
 if(flg_pro==TRUE){
  print(.message_pro)
 }
 
 return(res)
 
}

markov_model=function(Data, var_path, var_conv, var_value=NULL, var_null=NULL, order=1, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", ncore=1, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, verbose=TRUE, flg_pro=TRUE){
 
 
 if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
  print("Data must be a data.frame or a data.table")
 } 
 
 if(is.character(var_path)){
	if(!var_path%in%names(Data)){
	 print("var_path must be a column of Data")
	}
 }else{
  print("var_path must be a string")
 }
 
 if(is.character(var_conv)){
	if(!var_conv%in%names(Data)){
	 print("var_conv must be a column of Data")
	}
 }else{
  print("var_conv must be a string")
 }
 
 if(!is.null(var_value)){
  if(!var_value%in%names(Data)){
   print("var_value must be a column of Data")
  }
 }
 
 if(!is.null(var_null)){
  if(!var_null%in%names(Data)){
   print("var_null must be a column of Data")
  }
 }
 
 if(order<1){stop("order must be >= 1")}
 if(nsim_start<1){stop("nsim_start must be >= 1")}
 if(!is.null(max_step)){if(max_step<1){stop("max_step must be >= 1")}}
 if(!out_more%in%c(0,1)){stop("out_more must be FALSE or TRUE")}
 if(length(sep)>1){stop("sep must have length 1")}
 if(ncore<1){stop("ncore must be >= 1")}
 if(nfold<1){stop("nfold must be >= 1")}
 if(seed<0){stop("seed must be >= 0")}
 if(conv_par<0){stop("conv_par must be > 0")}
 if(rate_step_sim<0){stop("rate_step_sim must be > 0")}
 if(!verbose%in%c(0,1)){stop("verbose must be FALSE or TRUE")}
 
 if(nrow(Data[which(Data[var_conv]!=0),])==0){stop("Data must have at least one converting path")}

 if(is.null(var_value)){var_value="0"}
 if(is.null(var_null)){var_null="0"}
 if(is.null(max_step)){max_step=0}
 if(!is.null(seed)){set.seed(seed)}
 
 res=.Call("markov_model_cpp", Data, var_path, var_conv, var_value, var_null, order, nsim_start, max_step, out_more, sep, ncore, nfold, seed, conv_par, rate_step_sim,verbose)
 
 if(flg_pro==TRUE){
  print(.message_pro)
 }
 
 if(out_more==FALSE){
  return(as.data.frame(res)) 
 }else{
  return(list(result=as.data.frame(res$result),transition_matrix=as.data.frame(res$transition_matrix),removal_effects=as.data.frame(res$removal_effects)))
 }

}
 
 
transition_matrix=function(Data, var_path, var_conv, var_null, order=1, sep=">", flg_equal=TRUE, flg_pro=TRUE){
 
 if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
  print("Data must be a data.frame or a data.table")
 } 
 
 if(is.character(var_path)){
	if(!var_path%in%names(Data)){
	 print("var_path must be a column of Data")
	}
 }else{
  print("var_path must be a string")
 }
 if(is.character(var_conv)){
	if(!var_conv%in%names(Data)){
	 print("var_conv must be a column of Data")
	}
 }else{
  print("var_conv must be a string")
 }
  
 if(!is.null(var_null)){
  if(!var_null%in%names(Data)){
   print("var_null must be a column of Data")
  }
 }
 
 if(order<1){stop("order must be >= 1")}
 if(length(sep)>1){stop("sep must have length 1")}
 if(!flg_equal%in%c(0,1)){stop("flg_equal must be FALSE or TRUE")}

 if(is.null(var_null)){var_null="0"}

 res=.Call("transition_matrix_cpp", Data, var_path, var_conv, var_null, order, sep, flg_equal)
 
 if(flg_pro==TRUE){
  print(.message_pro)
 }
 
 return(list(channels=data.frame(id=1:length(res$channels),channel_name=res$channels),transition_matrix=as.data.frame(res$transition_matrix)))
  
}
 
 
auto_markov_model=function(Data, var_path, var_conv, var_null, var_value=NULL, max_order=10, roc_npt=100, plot=FALSE, nsim_start=1e5, max_step=NULL, out_more=FALSE, sep=">", ncore=1, nfold=10, seed=0, conv_par=0.05, rate_step_sim=1.5, verbose=TRUE, flg_pro=TRUE){
 
 if(!("data.frame"%in%class(Data)|"data.table"%in%class(Data))){
  print("Data must be a data.frame or a data.table")
 } 
 
 if(is.character(var_path)){
	if(!var_path%in%names(Data)){
	 print("var_path must be a column of Data")
	}
 }else{
  print("var_path must be a string")
 }
 if(is.character(var_conv)){
	if(!var_conv%in%names(Data)){
	 print("var_conv must be a column of Data")
	}
 }else{
  print("var_conv must be a string")
 }
 
 if(!is.null(var_value)){
  if(!var_value%in%names(Data)){
   print("var_value must be a column of Data")
  }
 }
 
 if(!is.null(var_null)){
  if(!var_null%in%names(Data)){
   print("var_null must be a column of Data")
  }
 }
 
 if(max_order<1){stop("max_order must be >= 1")}
 if(roc_npt<10){stop("roc_npt must be >= 10")}
 if(!plot%in%c(0,1)){stop("plot must be FALSE or TRUE")}
 if(nsim_start<1){stop("nsim_start must be >= 1")}
 if(!is.null(max_step)){if(max_step<1){stop("max_step must be >= 1")}}
 if(!out_more%in%c(0,1)){stop("out_more must be FALSE or TRUE")}
 if(length(sep)>1){stop("sep must have length 1")}
 if(ncore<1){stop("ncore must be >= 1")}
 if(nfold<1){stop("nfold must be >= 1")}
 if(seed<0){stop("seed must be >= 0")}
 if(conv_par<0){stop("conv_par must be > 0")}
 if(rate_step_sim<0){stop("rate_step_sim must be > 0")}
 if(!verbose%in%c(0,1)){stop("verbose must be FALSE or TRUE")}
 
 order=choose_order(Data, var_path, var_conv, var_null, max_order=max_order, sep=sep, ncore=ncore, roc_npt=roc_npt, plot=plot, flg_pro=FALSE)
 order=order[['suggested_order']]
 
 res=markov_model(Data, var_path, var_conv, var_value=var_value, var_null=var_null, order=order, nsim_start=nsim_start, max_step=max_step, out_more=out_more, sep=sep, ncore=ncore, nfold=nfold, seed=seed, conv_par=conv_par, rate_step_sim=rate_step_sim, verbose=verbose, flg_pro=FALSE)
 
 if(flg_pro==TRUE){
  print(.message_pro)
 }

 if(out_more==FALSE){
  return(as.data.frame(res)) 
 }else{
  return(list(result=as.data.frame(res$result),transition_matrix=as.data.frame(res$transition_matrix),removal_effects=as.data.frame(res$removal_effects)))
 }

}


.request_token_channelattributionpro = function(
    email,
    endpoint   = "https://app.channelattribution.io/genpkg/generate_token.php",
    timeout    = 10,
    verify_ssl = TRUE
  ) {

    # Basic email validation (simple but effective for most cases)
    is_valid_email = function(x) {
      is.character(x) && length(x) == 1L && nzchar(x) &&
        grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", x)
    }
    if (!is_valid_email(email)) {
      stop("Please enter a valid, non-empty email address (e.g., john.black@company.com).", call. = FALSE)
    }

    # Build a curl handle with common options
    make_handle = function(method = c("POST", "GET")) {
      method = match.arg(method)
      h = curl::new_handle()
      curl::handle_setheaders(h,
        "User-Agent" = "capro-token-client/1.0"
      )
      curl::handle_setopt(
        h,
        followlocation = TRUE,
        # total timeout (seconds)
        timeout = as.numeric(timeout)
      )
      # SSL verification knobs
      curl::handle_setopt(
        h,
        ssl_verifypeer = isTRUE(verify_ssl),
        ssl_verifyhost = if (isTRUE(verify_ssl)) 2L else 0L
      )
      if (identical(method, "POST")) {
        curl::handle_setform(h, email = email)
      }
      h
    }
  
    # Helper to perform request and capture status + body safely
    safe_fetch = function(url, handle, query = NULL) {
      full_url = if (is.null(query)) url else {
        paste0(url, if (grepl("\\?", url, fixed = TRUE)) "&" else "?", curl::curl_escape(names(query)), "=", curl::curl_escape(query))
      }
      out = tryCatch({
        res = curl::curl_fetch_memory(full_url, handle)
        list(
          ok     = TRUE,
          status = res$status_code,
          body   = rawToChar(res$content %||% raw())
        )
      }, error = function(e) {
        # Distinguish likely TLS/timeout vs generic transport
        msg = conditionMessage(e)
        if (grepl("SSL|TLS|certificate|timeout|timed out", msg, ignore.case = TRUE)) {
          stop(paste0("network_or_ssl_error: ", msg), call. = FALSE)
        } else {
          stop(paste0("request_error: ", msg), call. = FALSE)
        }
      })
      out$body = trimws(out$body %||% "")
      out
    }
  
    `%||%` = function(x, y) if (is.null(x)) y else x

    # 1) Try POST
    post_h = make_handle("POST")
    post_res = safe_fetch(endpoint, post_h)
  
    # 2) If 405/403 and body hints at method, retry with GET
    if (post_res$status %in% c(405L, 403L) &&
        nzchar(post_res$body) &&
        grepl("method", post_res$body, ignore.case = TRUE)) {
  
      get_h = make_handle("GET")
      get_res = safe_fetch(endpoint, get_h, query = c(email = email))
      return(get_res$body)
    }
  
    # Return body regardless of HTTP status
    post_res$body
  }


install_pro = function() {

  # Base-R secret prompt (no dependencies, works on Unix + Windows)
  .read_secret <- function(prompt = "Enter value: ") {
   trimws(readline(prompt))
  }

  if (!requireNamespace("curl", quietly = TRUE) || !requireNamespace("jsonlite", quietly = TRUE)) {
    stop("This function requires 'curl' and 'jsonlite'. Please install them from CRAN.", call. = FALSE)
  }
  
  # Prompt: token or email
  msg = "Enter your ChannelAttributionPro token. If you don't have one, enter your work/university email to request it: "
  token = .read_secret(msg)

  token = trimws(token)
  
  # If it looks like an email, trigger the token request and stop
  if (grepl("@", token, fixed = TRUE)) {
    email=token
    # simple email sanity check
    if (!grepl("^[^@\\s]+@[^@\\s]+\\.[^@\\s]+$", token)) {
      stop("Please enter a valid email address or a token.", call. = FALSE)
    }
    message("Sending a token...")
    res_token=.request_token_channelattributionpro(email=email)
    message("*** We email the token to eligible work or university addresses - check your inbox and Spam/Junk; if you don't receive it, try a different work/university email, and if it still doesn't arrive, contact info@channelattribution.io.")
    return(invisible(NULL))
  }
  
  if (!nzchar(token)) stop("A non-empty token or email is required.", call. = FALSE)

  `%||%` = function(x, y) if (is.null(x)) y else x

  urlencode_form = function(lst) {
    # Build application/x-www-form-urlencoded body
    if (length(lst) == 0) return("")
    keys = names(lst)
    if (is.null(keys)) stop("Form list must be named", call. = FALSE)
    kv = character(length(lst))
    for (i in seq_along(lst)) {
      k = curl::curl_escape(keys[i])
      v = curl::curl_escape(as.character(lst[[i]]))
      kv[i] = paste0(k, "=", v)
    }
    paste(kv, collapse = "&")
  }

  os_release_value = function(keys = c("ID_LIKE", "ID")) {
    path = "/etc/os-release"
    if (!file.exists(path)) return("")
    txt = try(readLines(path, warn = FALSE, encoding = "UTF-8"), silent = TRUE)
    if (inherits(txt, "try-error")) return("")
    vals = list()
    for (ln in txt) {
      if (!nzchar(ln) || grepl("^\\s*#", ln)) next
      kv = strsplit(ln, "=", fixed = TRUE)[[1]]
      if (length(kv) == 2L) {
        k = trimws(kv[1]); v = gsub('^"|"$', "", trimws(kv[2])); vals[[k]] = v
      }
    }
    for (k in keys) if (!is.null(vals[[k]])) return(vals[[k]])
    ""
  }

  # ---- POST helpers now send x-www-form-urlencoded (not multipart) ----
  http_post_form = function(url, form, timeout = 300) {
    h = curl::new_handle()
    curl::handle_setheaders(h,
      "User-Agent" = "capro-r-installer/1.1",
      "Accept" = "application/json,text/html;q=0.8,*/*;q=0.5",
      "Content-Type" = "application/x-www-form-urlencoded"
    )
    curl::handle_setopt(h, timeout = as.numeric(timeout), followlocation = TRUE)
    curl::handle_setopt(h, ssl_verifypeer = TRUE, ssl_verifyhost = 2L)
    body = urlencode_form(form)
    curl::handle_setopt(h, postfields = body)
    out = try(curl::curl_fetch_memory(url, handle = h), silent = TRUE)
    if (inherits(out, "try-error")) {
      msg = conditionMessage(attr(out, "condition"))
      body = if (grepl("SSL|TLS|certificate|timeout|timed out", msg, ignore.case = TRUE)) {
        paste0("network_or_ssl_error: ", msg)
      } else paste0("request_error: ", msg)
      return(list(status = 0L, body = body, headers = list()))
    }
    list(
      status  = as.integer(out$status_code),
      body    = rawToChar(out$content %||% raw()),
      headers = out$headers %||% list()
    )
  }

  http_get_text = function(url, timeout = 60) {
    h = curl::new_handle()
    curl::handle_setheaders(h, "User-Agent" = "capro-r-installer/1.1")
    curl::handle_setopt(h, timeout = as.numeric(timeout), followlocation = TRUE)
    curl::handle_setopt(h, ssl_verifypeer = TRUE, ssl_verifyhost = 2L)
    out = try(curl::curl_fetch_memory(url, handle = h), silent = TRUE)
    if (inherits(out, "try-error")) return(NULL)
    rawToChar(out$content %||% raw())
  }

  extract_href_links = function(html_text) {
    if (is.null(html_text) || !nzchar(html_text)) return(character())
    m = gregexpr("<a[^>]+href\\s*=\\s*\"([^\"]+)\"", html_text, ignore.case = TRUE, perl = TRUE)
    hits = regmatches(html_text, m)[[1]]
    if (!length(hits)) return(character())
    sub("\".*$", "", sub("^.*href\\s*=\\s*\"", "", hits, perl = TRUE), perl = TRUE)
  }

  resolve_pkg_url = function(pkg_value) {
    if (is.null(pkg_value)) return(NULL)
    pkg_value = trimws(as.character(pkg_value))
    if (!nzchar(pkg_value)) return(NULL)
    if (grepl("\\.(zip|tgz|tar\\.gz)$", pkg_value, ignore.case = TRUE)) return(pkg_value)
    dir_url = paste0(sub("/+$", "", pkg_value), "/")
    listing = http_get_text(dir_url, timeout = 60)
    if (is.null(listing)) return(NULL)
    links = extract_href_links(listing)
    files = links[!grepl("/$", links) & !links %in% c("/", "../")]
    if (!length(files)) return(NULL)
    cand = files[grepl("\\.(zip|tgz|tar\\.gz)$", files, ignore.case = TRUE)]
    if (!length(cand)) cand = files
    cand = cand[order(cand)]
    paste0(dir_url, cand[length(cand)])
  }

  install_from_url = function(url, os) {
    type = if (identical(os, "windows") || identical(os, "macos")) "binary" else "source"
    message("Installing from: ", url)
    out = try(utils::install.packages(url, repos = NULL, type = type, quiet = FALSE), silent = TRUE)
    !inherits(out, "try-error")
  }

  system_info_list = function() {
    sys = Sys.info()
    comp = NA_character_
    try({
      v = suppressWarnings(system("R CMD config CC", intern = TRUE))
      if (length(v) && nzchar(v[1])) comp = v[1]
    }, silent = TRUE)
    if (!nzchar(comp)) {
      for (cmd in c("gcc --version", "clang --version")) {
        try({
          v = suppressWarnings(system(cmd, intern = TRUE, ignore.stderr = TRUE))
          if (length(v) && nzchar(v[1])) { comp = v[1]; break }
        }, silent = TRUE)
      }
    }
    list(
      os = sys[["sysname"]],
      release = sys[["release"]],
      machine = sys[["machine"]],
      r_version = paste0(R.version$major, ".", R.version$minor),
      compiler = if (nzchar(comp)) comp else "not found"
    )
  }

  system_info_json = function() {
    jsonlite::toJSON(system_info_list(), auto_unbox = TRUE, pretty = TRUE)
  }

  notify_package_request_once = function(token, action, info,
                                          endpoint = "https://app.channelattribution.io/genpkg/build_check_email.php",
                                          timeout = 10) {
    if (!nzchar(token)) return("missing_token_param")
    send = function(payload) {
      h = curl::new_handle()
      curl::handle_setheaders(h,
        "User-Agent" = "capro-build-check-r/1.1",
        "Accept" = "text/plain, */*",
        "Content-Type" = "application/x-www-form-urlencoded"
      )
      curl::handle_setopt(h, timeout = as.numeric(timeout), followlocation = TRUE)
      curl::handle_setopt(h, ssl_verifypeer = TRUE, ssl_verifyhost = 2L)
      body = urlencode_form(list(token = token, action = action, info = payload))
      curl::handle_setopt(h, postfields = body)
      out = try(curl::curl_fetch_memory(endpoint, handle = h), silent = TRUE)
      if (inherits(out, "try-error")) {
        msg = conditionMessage(attr(out, "condition"))
        if (grepl("SSL|TLS|certificate|timeout|timed out", msg, ignore.case = TRUE)) {
          return(paste0("network_or_ssl_error: ", msg))
        }
        return(paste0("request_error: ", msg))
      }
      trimws(rawToChar(out$content %||% raw()))
    }
    resp = send(info)
    if (grepl("network_or_ssl_error:|request_error:", resp)) {
      small = if (nchar(info, type = "bytes") > 1024L) substr(info, 1L, 1024L) else info
      resp2 = send(small)
      return(paste0(resp, " | retry: ", resp2))
    }
    resp
  }

  # ---------- final notifier ----------
  action = "ERROR"
  info_blob = ""
  notifier_response = NULL
  on.exit({
    tok = if (is.null(token)) Sys.getenv("CHPRO_TOKEN", "") else token
    info_to_send = info_blob
    if (nzchar(info_to_send) && nchar(info_to_send, type = "bytes") > 8192L) {
      info_to_send = paste0(substr(info_to_send, 1L, 8192L), "...(truncated)")
    }
    notifier_response <= try(notify_package_request_once(tok, action, info_to_send), silent = TRUE)
    msg = if (inherits(notifier_response, "try-error")) conditionMessage(attr(notifier_response, "condition")) else as.character(notifier_response)
    #message("[notifier] build_check_email.php response: ", msg)
    if (interactive()) try(utils::flush.console(), silent = TRUE)
  }, add = TRUE)

  # ---------- token ----------
  if (is.null(token) || !nzchar(token)) token = Sys.getenv("CHPRO_TOKEN", "")
  if (!nzchar(token)) {
    message("Missing token. Pass token=... or set CHPRO_TOKEN in the environment.")
    return(invisible(NULL))
  }

  # ---------- env detection ----------
  sysname = Sys.info()[["sysname"]]
  machine = tolower(Sys.info()[["machine"]] %||% R.version$arch %||% "")
  arch = if (grepl("x86_64|amd64", machine)) "amd64" else if (grepl("aarch64|arm64", machine)) "arm64" else "amd64"

  os = NA_character_
  os_vers = NA_character_
  if (identical(sysname, "Darwin")) {
    os = "macos"; os_vers = if (identical(arch, "amd64")) "13" else "15"
  } else if (identical(sysname, "Windows")) {
    os = "windows"; os_vers = "11"
  } else {
    id_like = os_release_value(c("ID_LIKE", "ID"))
    if (grepl("rhel|fedora|centos|rocky|almalinux", id_like, ignore.case = TRUE)) {
      os = "rhel"; os_vers = "8"
    } else {
      os = "ubuntu"; os_vers = "20"
    }
  }

  r_version_minor = paste0(R.version$major, ".", strsplit(R.version$minor, "\\.")[[1]][1])
  params = list(
    os        = os,
    os_vers   = os_vers,
    arch      = arch,
    lang      = "r",
    lang_vers = r_version_minor,
    replace   = "0",
    uctr      = "0",
    token     = token
  )
  params1 = list(
    os        = os,
    os_vers   = os_vers,
    arch      = arch,
    lang      = "r",
    lang_vers = r_version_minor,
    replace   = "0",
    uctr      = "0"
  )

  message(sprintf("Detected -> os=%s os_vers=%s arch=%s R=%s", os, os_vers, arch, r_version_minor))
  message("Building the package. Estimated time: 0-30 minutes. Please wait...")

  # ---------- call builder ----------
  builder_url = "https://app.channelattribution.io/genpkg/genpkg.php"
  res = http_post_form(builder_url, params, timeout = 300)
  status = res$status
  text   = res$body

  if (identical(status, 401L)) {
    message("*** Token non valid or expired. Write to info@channelattribution.io.")
    action = "ERROR"
    info_blob = jsonlite::toJSON(list(
      reason = "invalid_token",
      builder_status = status,
      system = system_info_list(),
      params = params1
    ), auto_unbox = TRUE, pretty = TRUE)
    return(invisible(NULL))
  }

  data = NULL
  if (nzchar(text)) {
    tmp = try(jsonlite::fromJSON(text, simplifyVector = TRUE), silent = TRUE)
    if (!inherits(tmp, "try-error")) data = tmp
  }

  if (is.list(data)) {
    err  = tolower(as.character(data$error %||% ""))
    stat = tolower(as.character(data$status %||% ""))
    if (grepl("invalid token", err) || (stat %in% c("fail", "error") && grepl("token", err))) {
      message("Token non valid or expired. Write to info@channelattribution.io.")
      action = "ERROR"
      info_blob = jsonlite::toJSON(list(
        reason = "invalid_token_in_body",
        builder_status = status,
        body = substr(text, 1L, 800L),
        system = system_info_list(),
        params = params1
      ), auto_unbox = TRUE, pretty = TRUE)
      return(invisible(NULL))
    }
  }

  pkg_url = NULL
  ok_path = TRUE
  if (status %in% c(200L, 409L) && is.list(data) && !is.null(data$pkg)) {
    pkg_url = resolve_pkg_url(as.character(data$pkg)[1])
    if (is.null(pkg_url)) ok_path = FALSE
  } else {
    ok_path = FALSE
    message("Installation failed. Builder response:")
    message(sprintf("  HTTP %s", status))
    message(sprintf("  Body[0:800]: %s", substr(text %||% "", 1L, 800L)))
    message("Send the following information:")
    cat(system_info_json(), "\n\n")
    message("to info@channelattribution.io.")
    action = "ERROR"
    info_blob = jsonlite::toJSON(list(
      result = "builder_unexpected_response",
      builder_status = status,
      body = substr(text %||% "", 1L, 800L),
      system = system_info_list(),
      params = params1
    ), auto_unbox = TRUE, pretty = TRUE)
  }

  if (ok_path && nzchar(pkg_url)) {
    ok = isTRUE(install_from_url(pkg_url, os = os))
    if (ok) {
      message("*** Package installed. Restart the session and try: library(ChannelAttributionPro)")
      action = "SUCCESS"
      info_blob = jsonlite::toJSON(list(
        result = "installed",
        package = pkg_url,
        system = system_info_list(),
        params = params1
      ), auto_unbox = TRUE, pretty = TRUE)
      return(invisible(NULL))
    } else {
      message("Installation failed. Send the following information:")
      cat(system_info_json(), "\n\n")
      message("to info@channelattribution.io.")
      action = "ERROR"
      info_blob = jsonlite::toJSON(list(
        result = "install_failed",
        package = pkg_url,
        system = system_info_list(),
        params = params1
      ), auto_unbox = TRUE, pretty = TRUE)
      return(invisible(NULL))
    }
  }

  invisible(NULL)
}

