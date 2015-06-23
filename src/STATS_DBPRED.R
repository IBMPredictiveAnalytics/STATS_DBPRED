#Licensed Materials - Property of IBM
#IBM SPSS Products: Statistics General
#(c) Copyright IBM Corp. 2014
#US Government Users Restricted Rights - Use, duplication or disclosure 
#restricted by GSA ADP Schedule Contract with IBM Corp.

# Author: JKP, IBM SPSS
# Version = 1.0.1

# history
# 25-Jun-2013 - original version


helptext="STATS DBPRED
DSNAME=dataset name ID=variable name
WORKSPACE = file name
USERETAINEDWS = YES or NO*
/HELP

This procedure performs density-based clustering based
on results from the companion procedure STATS DBSCAN.

Example:
STATS DBPRED DSNAME=predicted WORKSPACE='c:/temp/clus.rdata'
ID = idvar.

The names of the variables to be used to predict the clusters
are taken from the saved or retained clustering workspace.
They must exist in the current active dataset.  Cases with
any missing data will not appear in the output dataset.

DSNAME specifies the name for the output dataset containing
the predicted clusters.  The dataset name must not already
be in use.

ID is the optional name of a variable in the active dataset to use as 
the ID in the output dataset.

Either a saved workspace file or a retained workspace containing
the output from the STATS DBSCAN must be specified.  Use
WORKSPACE to specify a file or specify USERETAINEDWS=YES
to use the most recent workspace retained in the clustering 
command.

STATS DBPRED /HELP prints this help and does nothing else.
"

spssdbscanpred <- function(dsname, id=NULL, useretainedws=FALSE, workspace=NULL) {
    # main routine

    setuplocalization("STATS_DBPRED")

    # A warnings proc name is associated with the regular output
    # (and the same omsid), because warnings/errors may appear in
    # a separate procedure block following the regular output
    procname=gtxt("Density-Based Clustering Prediction")
    omsid="STATSDBPRED"
    warns = Warn(procname=gtxt("Density-Based Clustering Prediction: Warnings"),
        omsid=omsid)

    tryCatch(library(fpc), error=function(e){
        warns$warn(gtxtf("The R %s package is required but could not be loaded.","fpc"),
            dostop=TRUE)
    }
    )

    # make sure active dataset has a name
    # and name is not in use
    alldsspecs = c(dsname)
    if (!is.null(alldsspecs)) {
        alldatasets = spssdata.GetDataSetList()
        if ("*" %in% alldatasets) {
            warns$warn(gtxt("The active dataset must have a name"),
                dostop=TRUE)
        }
        if (length(intersect(alldsspecs, alldatasets) > 0)) {
            warns$warn(gtxt("The output dataset name is already in use"),
                dostop=TRUE)
        }
    }
    if (!xor(useretainedws, !is.null(workspace))) {
        warns$warn(gtxt("Either a workspace file or the use-retained option but not both must be specified"), dostop=TRUE)
    }
        
    # If model is from a file, load it now
    if (!is.null(workspace)) {
        load(file=workspace)
    }
    # Is it a valid workspace?
    if (!(exists("stats_dbscan_res") && exists("stats_dbscan_details") 
        && exists("stats_dbscan_dta"))) {
        warns$warn(gtxt("One or more required objects are missing from the specified workspace"),
            dostop=TRUE)
    }
    if (is.null(stats_dbscan_res$isseed)) {
        warns$warn(gtxt("The specified workspace does not contain seeds.  It cannot be used for prediction."),
            dostop=TRUE)
    }
    # validate input dataset and fetch variables for prediction and id
    # Variable types are not checked.  Mismatch will cause prediction failure.
    
    alldata = c(stats_dbscan_details[["variables"]], id)
    vardict = spssdictionary.GetDictionaryFromSPSS()
    existingnames = vardict["varName",]
    missingvars = setdiff(stats_dbscan_details[["variables"]], existingnames)
    if (length(missingvars) > 0) {
        warns$warn(gtxtf("The following required variables are not in the active dataset: %s",
            paste(missingvars, collapse=" ")), dostop=TRUE)
    }
    dta = spssdata.GetDataFromSPSS(alldata, missingValueToNA=TRUE)
    nvars = length(stats_dbscan_details[["variables"]])
    
    # dbscan does not support factors or missing data.  Discard 
    # missing cases listwise (including missing id values).
    allcasecount = nrow(dta)
    dta = dta[complete.cases(dta),]
    completecasecount = nrow(dta)
    if (!is.null(id)) {
        iddata= dta[nvars+1]
        dta = dta[-(nvars+1)]
    } else {
        iddata = NULL
    }
    # Calculate predicted cluster.  Note that predict requires
    # the estimation data as well as the new data.
    predvalues = tryCatch(
        predict(
            stats_dbscan_res,
            stats_dbscan_dta,
            dta
        ),
        error = function(e) {
            warns$warn(e$message, dostop=TRUE)
            }
    )
        
    # print summary and create output dataset
    # 
    displayresults(warns, stats_dbscan_details, dsname, useretainedws, workspace,
        allcasecount, completecasecount)
    warns$display(inproc=TRUE)  # flushes warnings and ends procedure state
    createdataset(predvalues, iddata, stats_dbscan_details, dsname, vardict)

    # clean up workspace as needed
    rm(predvalues, dta)
    if (!useretainedws) {
        res <- tryCatch(rm(list=ls()), warning = function(e) {return(NULL)})
    }
}

displayresults = function(warns, details, dsname, useretainedws, workspace,
    allcasecount, completecasecount) {
    # Display summary of cluster parameters and prediction parameters
    
    StartProcedure(gtxt("Density-Based Clustering Prediction"), "STATSDBSCANPRED")

    # summary results
    scaption = gtxt("Computations done by R package fpc function dbscan by Christian Hennig")
    lbls = c(
        gtxt("Variables"), 
        gtxt("Method"),
        gtxt("Min Reachability Distance"),
        gtxt("Original Data Scaled"),
        gtxt("Minimum Cluster Size"),
        gtxt("Prediction Dataset"), 
        gtxt("Workspace File Used"),
        gtxt("Workspace Creation Date"),
        gtxt("Number of Prediction Cases"),
        gtxt("Number of Valid Prediction Cases")
    )
    vals = c(
        paste(details[["variables"]], collapse=" "),
        details["method"],
        details["rdist"],
        ifelse(details["scale"], gtxt("Yes"), gtxt("No")),
        details["minpts"],
        dsname,
        ifelse(useretainedws, gtxt("--Retained--"), workspace),
        details[["creationdate"]],
        allcasecount,
        completecasecount
    )

    # settings and result summary
    spsspivottable.Display(
        data.frame(cbind(vals), row.names=lbls), 
        title = gtxt("Settings and Results Summary"),
        collabels=c(gtxt("Summary")), 
        templateName="DBSCANPREDSUMMARY", 
        outline=gtxt("Summary"),
        caption = scaption
    )
}

createdataset = function(pred, iddata, details, dsname, vardict) {
    # Create classification dataset
    # Dataset name is known to be okay, and procedure state is ended
    # details contents are from the original clustering

    if (!is.null(iddata)) {  # was an id variable provided - copy properties
        pred = data.frame(iddata, pred)
        loc = match(details[["id"]], vardict["varName",])
        idtype = as.integer(vardict[["varType", loc]])
        idformat = vardict[["varFormat", loc]]
        idlabel = vardict[["varLabel", loc]]
    } else {
        pred = data.frame(id=1:nrow(pred), pred)
        idtype = 0
        idformat = "F10.0"
        idlabel = ""
    }

    dict = spssdictionary.CreateSPSSDictionary(
        c("ID", idlabel, idtype, idformat, "nominal"),
        c("PredCluster", gtxt("Predicted Cluster"), 0, "F8.0", "nominal")
    )
    spssdictionary.SetDictionaryToSPSS(dsname, dict)
    spssdata.SetDataToSPSS(dsname, pred)
    spssdictionary.EndDataStep()
}  
    
# localization initialization
setuplocalization = function(domain) {
    # find and bind translation file names
    # domain is the root name of the extension command .R file, e.g., "SPSSINC_BREUSCH_PAGAN"
    # This would be bound to root location/SPSSINC_BREUSCH_PAGAN/lang

    fpath = Find(file.exists, file.path(.libPaths(), paste(domain, ".R", sep="")))
    bindtextdomain(domain, file.path(dirname(fpath), domain, "lang"))
} 

# override for api to account for extra parameter in V19 and beyond
StartProcedure <- function(procname, omsid) {
    if (substr(spsspkg.GetSPSSVersion(),1, 2) >= 19) {
        spsspkg.StartProcedure(procname, omsid)
    }
    else {
        spsspkg.StartProcedure(omsid)
    }
}
gtxt <- function(...) {
    return(gettext(...,domain="STATS_DBPRED"))
}

gtxtf <- function(...) {
    return(gettextf(...,domain="STATS_DBPRED"))
}

Warn = function(procname, omsid) {
    # constructor (sort of) for message management
    lcl = list(
        procname=procname,
        omsid=omsid,
        msglist = list(),  # accumulate messages
        msgnum = 0
    )
    # This line is the key to this approach
    lcl = list2env(lcl) # makes this list into an environment

    lcl$warn = function(msg=NULL, dostop=FALSE, inproc=FALSE) {
        # Accumulate messages and, if dostop or no message, display all
        # messages and end procedure state
        # If dostop, issue a stop.

        if (!is.null(msg)) { # accumulate message
            assign("msgnum", lcl$msgnum + 1, envir=lcl)
            # There seems to be no way to update an object, only replace it
            m = lcl$msglist
            m[[lcl$msgnum]] = msg
            assign("msglist", m, envir=lcl)
        } 

    if (is.null(msg) || dostop) {
        lcl$display(inproc)  # display messages and end procedure state
        if (dostop) {
            stop(gtxt("End of procedure"), call.=FALSE)  # may result in dangling error text
        }
    }
}

    lcl$display = function(inproc=FALSE) {
        # display any accumulated messages as a warnings table or as prints
        # and end procedure state, if any

    if (lcl$msgnum == 0) {   # nothing to display
        if (inproc) {
            spsspkg.EndProcedure()
        }
    } else {
        if (!inproc) {
            procok =tryCatch({
                StartProcedure(lcl$procname, lcl$omsid)
                TRUE
                },
                error = function(e) {
                    FALSE
                }
            )
        }
        if (procok) {  # build and display a Warnings table if we can
            table = spss.BasePivotTable("Warnings ","Warnings") # do not translate this
            rowdim = BasePivotTable.Append(table,Dimension.Place.row, 
                gtxt("Message Number"), hideName = FALSE,hideLabels = FALSE)

    for (i in 1:lcl$msgnum) {
        rowcategory = spss.CellText.String(as.character(i))
        BasePivotTable.SetCategories(table,rowdim,rowcategory)
        BasePivotTable.SetCellValue(table,rowcategory, 
            spss.CellText.String(lcl$msglist[[i]]))
    }
    spsspkg.EndProcedure()   # implies display
} else { # can't produce a table
    for (i in 1:lcl$msgnum) {
        print(lcl$msglist[[i]])
    }
}
}
}
return(lcl)
}
Run<-function(args){
    
    cmdname = args[[1]]
    args <- args[[2]]
    oobj<-spsspkg.Syntax(templ=list(
        spsspkg.Template("DSNAME", subc="", ktype="varname", var="dsname"),
        spsspkg.Template("ID", subc="", ktype="existingvarlist", var="id"),
        spsspkg.Template("USERETAINEDWS", subc="", ktype="bool", var="useretainedws"),
        spsspkg.Template("WORKSPACE", subc="", ktype="literal", var="workspace")
    ))        
if ("HELP" %in% attr(args,"names"))
    #writeLines(helptext)
    helper(cmdname)
else
    res <- spsspkg.processcmd(oobj,args,"spssdbscanpred")
}

helper = function(cmdname) {
    # find the html help file and display in the default browser
    # cmdname may have blanks that need to be converted to _ to match the file
    
    fn = gsub(" ", "_", cmdname, fixed=TRUE)
    thefile = Find(file.exists, file.path(.libPaths(), fn, "markdown.html"))
    if (is.null(thefile)) {
        print("Help file not found")
    } else {
        browseURL(paste("file://", thefile, sep=""))
    }
}
if (exists("spsspkg.helper")) {
assign("helper", spsspkg.helper)
}

