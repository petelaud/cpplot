
# Create contour plot for a single cell of the panel plot
CPcontour <- function(plotdata,
                      alpha,
                      par3 = 0,
                      nums,
                      xlim = c(0,1),
                      ylim = xlim,
                      methlab = "MN",
                      avg = F,
                      lside = F,
                      lines = F,
                      lines1 = F,
                      CIlen = F,
                      locind = F,
                      crude = F,
                      res.factor = 6,
                      colour = T,
                      textsize = 1) {

  if (colour == T) {
  	col1 = colorRampPalette(c("BLACK","RED"))(5)[4]
  	col2 = "DARKORANGE"
  	col3 = "#FFBF00"
  	col4 = colorRampPalette(c("WHITE","YELLOW"))(5)[4]
  	contcol = "BLACK"
  	contwd = 0.75
  } else {
  	col1="BLACK"
  	col2="GRAY33"
  	col3="GRAY66"
  	col4="WHITE"
  	contcol="GRAY25"
  	contwd=0.6
  }

	cols <- switch(as.character(alpha),
		"0.2" = c(rep("#000000",120),
		        rep(colorRampPalette(c("BLACK",col1))(10)[-1],each=4),
		        rep(c(col2,col3),each=4),
		        rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(9),each=4)),
		"0.1" = c(rep("#000000",160),
		        rep(colorRampPalette(c("BLACK",col1))(10)[-1],each=2),
		        col2,col2,col3,col3,
		        rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(9),each=2)),
		"0.05" = (c(rep("#000000",180),
		          colorRampPalette(c("BLACK",col1))(10)[-1],
		          col2,col3,colorRampPalette(c(col4,"WHITE"),space="Lab")(9))),
		"0.01" = c(rep("#000000",1960),
		         rep(colorRampPalette(c("BLACK",col1))(10)[-1],each=2),
		         rep(c(col2,col3),each=2),
		         rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(9),each=2))
	)
	oscols <- switch(as.character(alpha),
	  "0.5" = rev(c(rep("#000000",160),
	                rep(colorRampPalette(c("BLACK",col1))(10)[-1],each=2),
	                col2,col2,col3,col3,
	                rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(9),each=2))),
	  "0.2" = rev(c(rep("#000000",160),
		              rep(colorRampPalette(c("BLACK",col1))(10)[-1],each=2),
		              col2,col2,col3,col3,
		              rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(9),each=2))),
		"0.1" = rev(c(rep("#000000",180),
		              rep(colorRampPalette(c("BLACK",col1))(9)[-1],each=1),
                  col2,col2,col3,col3,
		              rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(8),each=1))),
		"0.05" = rev(c(rep("#000000",380),
		               rep(colorRampPalette(c("BLACK",col1))(9)[-1],each=1),
		               col2,col2,col3,col3,
		               rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(8),each=1))),
		"0.01" = rev(c(rep("#000000",1980),
		               rep(colorRampPalette(c("BLACK",col1))(9)[-1],each=1),
		               col2,col2,col3,col3,
		               rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(8),each=1)))
	)
	lencols1 <- c(colorRampPalette(c("WHITE",col4))(100),
	              colorRampPalette(c(col4,col1))(201))
	lencols <- c(colorRampPalette(c("WHITE","WHITE"))(100),
	             colorRampPalette(c("WHITE",col4))(100),
	             colorRampPalette(c(col1,"BLACK"))(101))
	loccols = (c(colorRampPalette(c("BLACK",col1))(5)[-1],
	              col2,col3,
	              colorRampPalette(c(col4,"WHITE"),space="Lab")(4)))

	conts <- switch(as.character(alpha),
		"0.05" = seq(0.005,0.995,0.005),
		"0.025" = seq(0,1,0.005),
		"0.1" = seq(0.005,0.995,0.005),
		"0.2" = seq(0.005,0.995,0.005),
		"0.01" = seq(0.001,0.999,0.001)
	)
	#optional alternative contours for 1-sided plots
	osconts <- switch(as.character(alpha),
		"0.05" = seq(0.0025,0.9975,0.0025),
		"0.025" = seq(0,1,0.0025),
		"0.1" = seq(0.0025,0.9975,0.0025),
		"0.2" = seq(0,1,0.005),
		"0.01" = seq(0.0005,0.9995,0.0005)
	)
	lenconts <- seq(-10,20,0.02)

	p1 <- as.numeric(dimnames(plotdata[["mastercp"]])[[1]])
	p2 <- as.numeric(dimnames(plotdata[["mastercp"]])[[2]])
	x <- p1[p1>=min(xlim) & p1<=max(xlim)]
	y <- p2[p2>=min(ylim) & p2<=max(ylim)]

	rawcpdata <- plotdata[["mastercp"]][paste(x),
	                                    paste(y),
	                                    paste(par3),
	                                    methlab,
	                                    paste(100*(1-alpha)),
	                                    "cp",
	                                    nums,]
	avecpdata <- ifelse(is.na(rawcpdata), NA, plotdata[["mastercp"]][paste(x),
	                                    paste(y),
	                                    paste(par3),
	                                    methlab,
	                                    paste(100*(1-alpha)),
	                                    "avecp",
	                                    nums,])
	if(avg) {
		if(lside == TRUE) cpdata <- plotdata[["mastercp"]][paste(x),
		                                          paste(y),
		                                          paste(par3),
		                                          methlab,
		                                          paste(100*(1-alpha)),
		                                          "averncp",
		                                          nums,]
		else cpdata <- avecpdata
	} else {
	  lirev <- 1 - plotdata[["mastercp"]][paste(x),
	                               paste(y),
	                               paste(par3),
	                               methlab,
	                               paste(100*(1-alpha)),
	                               "locindex",
	                               nums,]
	  ncpval <- 1 - plotdata[["mastercp"]][paste(x),
	                               paste(y),
	                               paste(par3),
	                               methlab,
	                               paste(100*(1-alpha)),
	                               "cp",
	                               nums,]

	  if(lside == TRUE) cpdata <- lirev * ncpval
	  else cpdata <- rawcpdata
	}
	if(CIlen) {
		if(avg) cpdata <- plotdata[["mastercp"]][paste(x),
		                                paste(y),
		                                methlab,
		                                paste(100*(1-alpha)),
		                                "len",
		                                nums,] - plotdata[["mastercp"]][,,"MN",
		                                                    paste(100*(1-alpha)),
		                                                    "len",
		                                                    nums,]
		else cpdata <- plotdata[["mastercp"]][paste(x),
		                             paste(y),
		                             methlab,
		                             paste(100*(1-alpha)),
		                             "len",
		                             nums,]
		cpdata[cpdata<(-10)] <- (-10)
		cpdata[cpdata>20] <- 20
	}
	if(locind) {
	  # Rename this
	  avecpdata <- plotdata[["mastercp"]][paste(x),
	                                           paste(y),
	                                           paste(par3),
	                                           methlab,
	                                           paste(100*(1-alpha)),
	                                           "avelocindex",
	                                           nums,]
	  if(avg) cpdata <- plotdata[["mastercp"]][paste(x),
	                                           paste(y),
	                                           paste(par3),
	                                           methlab,
	                                           paste(100*(1-alpha)),
	                                           "avelocindex",
	                                           nums,]
	  else cpdata <- plotdata[["mastercp"]][paste(x),
	                                        paste(y),
	                                        paste(par3),
	                                        methlab,
	                                        paste(100*(1-alpha)),
	                                        "locindex",
	                                        nums,]
	}

	if(!CIlen) cpdata[cpdata>1] <- 1
	if(alpha > 0.05) {
	  ccol <- ceiling(200*cpdata)  #colours for distinguishing areas above/below nominal level on the 3d surface plot
	} else if (alpha == 0.01) {
	  ccol <- ceiling(2000*cpdata)
	} else if (alpha == 0.025 | alpha == 0.05) {
	  ccol <- ceiling(200*cpdata*ifelse(lside,2,1))
	}
	if(CIlen) ccol <- ceiling(20*cpdata)
	if(locind) ccol <- ceiling(10*cpdata)

	if(!CIlen) ccol[ccol<1] <- 1

	z <- cpdata #[p1>=min(xlim) & p1<=max(xlim),p2>=min(ylim) & p2<=max(ylim)]
	zc <- avecpdata

	shift <- 0
	palette <- switch(as.character(lside), "TRUE" = oscols, "FALSE" = cols)
	if (CIlen) {
		if(avg) {
		  palette <- lencols
		  shift=200
		} else palette <- lencols1
	}
	if (locind) palette <- loccols
#	par(mar=res.factor*(c(2,1,2,0)+0.1))
	plot(x, y, type = 'n', col.axis="white", xlim=xlim,
	     ylim=ylim, 	      xaxs="i",
	     yaxs="i",)
	abline(h = seq(0, 1, 0.02), v = seq(0, 1, 0.025), col="black")
	image(x,
	      y,
	      ccol,
	      add = TRUE,
	      col = palette[min(ccol+shift,na.rm=T):max(ccol+shift,na.rm=T)],
	      xlab="",
	      ylab="",
	      cex.lab=res.factor*0.8*textsize,
	      col.axis="white",
	      xaxs="i",
	      yaxs="i",
	      tcl=0,
	      xlim=xlim,
	      ylim=ylim
	      )
	if (avg | CIlen) {
	if (locind == FALSE) {
	    contour(x,
	          y,
	          zc,
	          levels=switch(as.character(CIlen),
	                        "TRUE"=lenconts,
	                        "FALSE"=switch(as.character(lside),
	                                       "TRUE"=osconts,
	                                       "FALSE"=conts)),
	          add=T,
	          labcex=0.5*textsize*res.factor,
	          lwd=contwd*textsize*res.factor,
	          col=contcol,
	          vfont=c("sans serif","bold")
	          )
	}
}
	axis(side = 1,
	     at = seq(0,xlim[2],xlim[2]*0.1),
	     labels = c("0",rep("",9),paste(xlim[2])),
	     padj = 1,
	     lwd = 0.5*res.factor,
	     tcl = -1,
	     cex.axis = (2/3)*textsize*res.factor)
	mtext(side = 1,
	      text = bquote(italic(p)[1]),
	      cex = res.factor*0.75,
	      line = 1*res.factor)
	axis(side = 2,
	     at = seq(0,ylim[2],0.1*ylim[2]),
	     labels = c("0",rep("",9),paste(ylim[2])),
	     padj = 0,
	     lwd = 0.5*res.factor,
	     tcl = -1,
	     cex.axis = (2/3)*textsize*res.factor)
	mtext(side = 2,
	      text = bquote(italic(p)[2]),
	      cex = res.factor*0.75*textsize,
	      line = res.factor/3)
	box(lwd = 0.6*textsize*res.factor)

	if(lines != FALSE){
		if (lside && lines %in% c("both","1")) {
		  lines(x=c(0.5,0.85), y=c(0.6,0.95), lwd = 0.5*textsize*res.factor)
		}	else if (lside==FALSE & lines %in% c("RR")) {
  		  abline(coef = c(0, 1/2), lwd = 0.5*textsize*res.factor)
  		  abline(coef = c(0, 1/3), lwd = 0.5*textsize*res.factor)
  		  abline(coef = c(0, 1/4.5), lwd = 0.5*textsize*res.factor)
  		  abline(coef = c(0, 1/6), lwd = 0.5*textsize*res.factor)
		}	else if (lside==FALSE & lines %in% c("RD")) {
		  abline(coef = c(0.1, 1), lwd = 0.5*textsize*res.factor)
		  abline(coef = c(0.2, 1), lwd = 0.5*textsize*res.factor)
		}
	}

	palette
}


###
# New version of plotpanel, showing unsmoothed cp and location index together
# Note: Contains some redundant code left over from previous version, which
# had an option to plot either two-sided or on-sided coverage, where
# now we show both in the same panel
# Legend code is repetitive but gets the job done
plotpanel <- function(plotdata,
                      alpha,
                      par3,
                      sel,
                      oneside = F,
                      plotlab,
                      linesx = F,
                      smoothed = TRUE,
                      CIlen = F,
                      res.factor = 6,
                      fmt = "png",
                      colour = T,
                      textsize = 1,
                      limits = c(0,1),
                      sided = "R"
) {

  g <- "gamma"
  nmeth <- length(sel)
  nums <- dimnames(plotdata[["summaries"]])[[5]]
  contrast <- dimnames(plotdata[["summaries"]])[[6]]
  summaries <- plotdata[["summaries"]][paste(par3),,paste(100*(1-alpha)),,nums,]
  longlab <- methods <- dimnames(plotdata[["mastercp"]])[[4]]
  # Update labels as preference for nomenclature evolved since starting the code
  # and it would take too long to go back and re-run everything
  longlab[longlab=="SCAS"] <- "SCASu"
  longlab[longlab=="SCAS-bc"] <- "SCAS"
  longlab[longlab=="MOVER-c5"] <- "MOVER-c (\U0263=0.5)"
  longlab[longlab=="MOVER-c25"] <- "MOVER-c (\U0263=0.25)"
  longlab[longlab=="MOVER-c125"] <- "MOVER-c (\U0263=0.125)"
  longlab[longlab=="SCAS-c5"] <- "SCAS-c (\U0263=0.5)"
  longlab[longlab=="SCAS-c25"] <- "SCAS-c (\U0263=0.25)"
  longlab[longlab=="SCAS-c125"] <- "SCAS-c (\U0263=0.125)"

  longlab[longlab=="SCASp"] <- "T-SCASp(N-1)"
  longlab[longlab=="SCASpu"] <- "T-SCASp"
  longlab[longlab=="mid-p"] <- "T-midp"
  longlab[longlab=="Jeffreys"] <- "T-Jeffreys"
  longlab[longlab=="Wilson"] <- "T-Wilson"
  longlab[longlab=="SCASp-c125"] <- "T-SCASp-c125"
  longlab[longlab=="SCASp-c25"] <- "T-SCASp-c25"
  longlab[longlab=="SCASp-c5"] <- "T-SCASp-c5"
  longlab[longlab=="midp-c25"] <- "T-midp-c25"
  longlab[longlab=="Jeffreys-c25"] <- "T-Jeffreys-c25"
  longlab[longlab=="Wilson-c"] <- "T-Wilson-c"
  longlab[longlab=="C-P"] <- "T-CP"

  names(longlab) <- names(methods) <- methods
  n.grid <- dim(plotdata[["mastercp"]])[1]
  if (colour == F) collab=".bw" else collab=""

  # Set up to adjust number 4 rows of plots
  rows <- ifelse(smoothed == TRUE, 4, 3)

  # Select plot output format depending on journal requirements
  if (fmt=="tiff")  {
    tiff(file = paste(outpath,"_",fmt,"/","summary",
                      ifelse(sided=="L","L",""),
                      plotlab, nums, "_", 100*(1-alpha), "_", format(par3, nsmall=2),
                      ifelse(oneside, "os", ""), collab, ".tiff",
                      sep=""
    ),
    width = (120 * nmeth + 60) * res.factor,
    height = 4*rows * 38 * res.factor,
    type="quartz"
    )
  } else if (fmt=="png") {
    png(file = paste(outpath,
                     ifelse(sided=="L","L",""),
                     plotlab, nums, "_", 100*(1-alpha), "_", format(par3, nsmall=2),
                     ifelse(oneside, "os", ""), collab, ".png",
                     sep=""
    ),
    width = (115 * nmeth + 60) * res.factor,
#    height = (4*rows) * 38 * res.factor,
    height = (4*rows+0) * 38 * res.factor,
    type = "quartz"
    )
  }
  if (FALSE) { # Left-over code not checked for this application
  #else
    if(fmt=="eps") {
    setEPS()
    postscript(file = paste(outpath,"_","png","/","summary",
                            ifelse(sided=="L","L",""),plotlab,100*(1-alpha),"_",
                            nums,ifelse(oneside,"os",""),collab,".eps",
                            sep=""
    ),
    width=60,
    height=33
    )
  } else {
    png(file = paste(outpath,"_",fmt,"/",plotlab,
                     format(100*(1-alpha),scientific=F),"/","summary_",
                     ifelse(sided=="L","L",""),
                     sub(",","_",nums),ifelse(oneside,"os",""),collab,".png",
                     sep=""
    ),
    width=600*res.factor,
    height=400*res.factor
    )
  }
  }

  if (smoothed == TRUE) {
    matr <- matrix(c(1:(rows * nmeth), rows * nmeth + c(1, 1:(rows-1))),
                   rows, (nmeth + 1), byrow = FALSE)
  } else matr <- matrix(c(1:(rows * nmeth), rows * nmeth + c(1:rows)),
                          rows, (nmeth + 1), byrow = FALSE)
  layout(matr, widths = c(rep(2, nmeth), 1), heights = rep(4, rows))
  par(oma = res.factor*c(0,3,ifelse(fmt=="xxx",10,7),0), pty='s')
  par(cex.main = res.factor*0.8*textsize, cex.axis=res.factor*0.8*textsize)
  par(mar = res.factor*(c(2,1,1,0.5)+0.1))

  for(i in methods[sel]) {
    par(mar = res.factor*(c(2,1,1,0.5)+0.1))
    methlab <- longlab[i]
    # Run the unsmoothed plot and capture the output palette
    palette <- CPcontour(plotdata=plotdata,
                         alpha=alpha,
                         par3=par3,
                         nums=nums,
                         methlab=i,
                         lside=oneside,
                         avg = FALSE,
                         lines=linesx,
                         locind = FALSE,
                         CIlen=CIlen,
                         res.factor=res.factor,
                         colour=colour,
                         textsize=textsize,
                         xlim=limits)
    mtext(side = 3,
          cex = res.factor*0.8*textsize,
          line = 7*res.factor,
          text = latex2exp::TeX(paste0("\\textbf{",methlab,"}")))
    mtext(side = 3,
          cex = res.factor*0.8*textsize,
          line = 0.5*res.factor,
          text = paste0(
      ifelse(oneside,
             paste0("\n",ifelse(sided=="R","maxRNCP=","maxLNCP="),
                    summaries[i,"maxLNCP"]),
             paste0("\n","meanCP=",summaries[i,"meanCP"],
                    "\n","minCP=",summaries[i,"minCP"],
                    "\n","pctCons=",summaries[i,"pctCons"],"%"
             )
      ),
      ifelse(any(sel=="SCcc"),
             paste0(("\n"),"pctCons=",(summaries[i,paste("pctCons",ifelse(oneside,".1side",""),sep="")]),"%"),
             paste0("\n","within ±",
                    (ifelse(oneside,
                            format(0.1*alpha/2,scientific=F),
                            format(0.1*alpha,scientific=F))
                    ),
                    "=",
                    (summaries[i,paste("pctnear",ifelse(oneside,".1side",""),sep="")]),
                    "%",
                    "\n","below ",
                    format(1-1.1*alpha,scientific=F),
                    "=",(summaries[i,"pctBad"]),
                    "%"
             )
      )))

  if (smoothed == TRUE) {
    # Run the smoothed plot
    CPcontour(plotdata=plotdata,
                         alpha=alpha,
                         par3=par3,
                         nums=nums,
                         methlab=i,
                         lside=oneside,
                         avg = TRUE,
                         lines=linesx,
                         locind = FALSE,
                         CIlen=CIlen,
                         res.factor=res.factor,
                         colour=colour,
                         textsize=textsize,
                         xlim=limits)
    mtext(side = 3,
          cex = res.factor*0.8*textsize,
          line = 0.5*res.factor,
          text = (paste0(#longlab[i],
            "MACP cons=",summaries[i,"pctAveCons"],"%",
            ifelse(any(sel=="SCcc"),
                               paste0(("\n"),"pctCons=",(summaries[i,paste("pctCons",ifelse(oneside,".1side",""),sep="")]),"%"),
                               paste0("\n","within ±",
                                      (ifelse(oneside,
                                              format(0.1*alpha/2,scientific=F),
                                              format(0.1*alpha,scientific=F))
                                      ),"=",
                                      (summaries[i,paste("pctAvenear",ifelse(oneside,".1side",""),sep="")]),
                                      "%"
                               )
                        )))
    )
  }
    # Run the location index plot
    palette2 <- CPcontour(plotdata = plotdata,
              alpha = alpha,
              par3=par3,
              nums = nums,
              methlab = i,
              lside = oneside,
              lines = linesx,
              locind = TRUE,
              res.factor = res.factor,
              colour = colour,
              textsize = textsize,
              xlim = limits)

    mtext(side=3,
          cex=res.factor*0.8*textsize,
          line=0.5*res.factor,
          text = (paste0(
                   paste0("\n","mean loc. index=",summaries[i,"meanlocindex"]),
            "\n","within ±0.1=",
            (summaries[i,"pctgoodloc"]),
            "%" #,
          ))
    )
    # Run the DNCP plot
    palettex <- CPcontour(plotdata = plotdata,
                          alpha = alpha,
                          par3 = par3,
                          nums = nums,
                          methlab = i,
                          lside = TRUE,
                          lines = linesx,
                          res.factor = res.factor,
                          colour = colour,
                          textsize = textsize,
                          xlim = limits)
    mtext(side=3,
          cex=res.factor*0.8*textsize,
          line=0.5*res.factor,
          text = (paste0(
            "\n","meanDNCP",
            "=",(summaries[i,"meanDNCP"]),
            "\n","DNCP above ",
            format(1.2*alpha/2,scientific=F),
            "=",(summaries[i,"pctBad.DNCP"]),
#            "%",
#            "\n","DNCP within ±",
#            format(0.2*alpha/2,scientific=F),
#            "=",(summaries[i,"pctnear.1side"]),
            "%"
          )))
  }

  par(pty="m", mar=res.factor*c(2.35, 1.5, ifelse(fmt=="xxx", 2, 1.35), 3)+0.1)

  # Add legend for colour palette

  if (alpha == 0.05) {
    if (oneside == T) {
      image(y = (1:20-0.5)/20,
            z = matrix(1:20,nrow=1),
            col = rev(palette)[381:400],
            axes = F,
            ylab = "")
      text(x=1.2,
           y=(c(190:200)-190)/10,
           labels=(10:0)/200,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    } else {
      image(y=(1:20-0.5)/20,
            z=matrix(1:20,nrow=1),
            col=palette[181:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10, 0.45, 0.55),
           labels = c((90:100)/100, 0.945, 0.955),
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    }
  } else if (alpha == 0.1) {
    if (oneside == T) {
      image(y=(1:20-0.5)/20,
            z=matrix(181:200,nrow=1),
            col=rev(palette)[181:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(90:100)-90)/10,
           labels=(10:0)/100,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    } else {
      image(y=(1:40-0.5)/40,
            z=matrix(161:200,nrow=1),
            col=palette[161:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10, 0.45, 0.55),
           labels=c((40:50)/50, 0.89, 0.99),
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    }
  } else if (alpha == 0.2) {
    if (oneside == T) {
      image(y=(1:40-0.5)/40,
            z=matrix(161:200,nrow=1),
            col=rev(palette)[161:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(80:100)-80)/20,
           labels=(20:0)/100,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    } else {
      image(y=(1:80-0.5)/80,
            z=matrix(121:200,nrow=1),
            col=palette[121:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(30:50)-30)/20,
           labels=(30:50)/50,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    }
  } else if (alpha == 0.01) {
    if (oneside == T) {
      image(y=(1:20-0.5)/20,
            z=matrix(1981:2000,nrow=1),
            col=rev(palette)[1981:2000],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(195:200)-195)/5,
           labels=(5:0)/200,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    } else {
      image(y=(1:40-0.5)/40,
            z=matrix(1961:2000,nrow=1),
            col=palette[1961:2000],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10, 0.45, 0.55),
           labels = c(((5*98):(5*100))/(5*100), 0.985, 0.995),
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
    }
  }
  box(lwd=0.6*textsize*res.factor)
  mtext(text=ifelse(oneside,
                    ifelse(sided=="R","RNCP","LNCP"),"CP"),
        side = 3,
        at = 0,
        adj = 0.5,
        cex=res.factor*0.6,
        line=2*res.factor/3)


  # Add legend for location index
  image(y = (1:20-0.5)/20,
        z = matrix(1:20,nrow=1),
        col = (palette2),
        axes = F,
        ylab = "")
  text(x=1.2,
       y=(0:10)/10,
       labels=(0:10)/10,
       xpd=T,
       cex=res.factor*textsize,
       pos=4)
  text(x=c(-2),
       y=c(0),
       labels=c("Too mesial"),
       adj=c(0),
       xpd=T,
       srt=90,
       cex=res.factor*textsize)
  text(x=c(-2),
       y=c(1),
       labels=c("Too distal"),
       adj=c(1),
       xpd=T,
       srt=90,
       cex=res.factor*textsize)

  box(lwd=0.6*textsize*res.factor)
  mtext(text="MNCP/NCP",
        side=3,
        at = 0,
        adj = 0.5,
        cex=res.factor*0.6,
        line=2*res.factor/3)

  # Add legend for DNCP
  if (alpha == 0.05) {
      image(y = (1:20-0.5)/20,
            z = matrix(1:20,nrow=1),
            col = rev(palettex)[381:400],
            axes = F,
            ylab = "")
      text(x=1.2,
           y=(c(190:200)-190)/10,
           labels=(10:0)/200,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
  } else if (alpha == 0.1) {
      image(y=(1:20-0.5)/20,
            z=matrix(181:200,nrow=1),
            col=rev(palettex)[181:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(90:100)-90)/10,
           labels=(10:0)/100,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)

  } else if (alpha == 0.2) {
      image(y=(1:40-0.5)/40,
            z=matrix(161:200,nrow=1),
            col=rev(palettex)[161:200],
            axes=F,
            ylab="")
      text(x=1.2,
           y=(c(80:100)-80)/20,
           labels=(20:0)/100,
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
  } else if (alpha == 0.01) {
    image(y=(1:20-0.5)/20,
          z=matrix(1981:2000,nrow=1),
          col=rev(palettex)[1981:2000],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10),
           labels=(10:0)/(1000),
           xpd=T,
           cex=res.factor*textsize,
           pos=4)
  }
  box(lwd=0.6*textsize*res.factor)
  mtext(text="DNCP",
        side=3,
        at = 0,
        adj = 0.5,
        cex=res.factor*0.6,
        line=2*res.factor/3)

  if(fmt=="xxx") {
    mtext(paste(ifelse(oneside,
                       "Left-side ",""
    ),
    "Coverage properties of ",100*(1-alpha),"%"," CIs for ","p1-p2",
    " with (n1,n2)=(",nums,")",
    "\n",
    ifelse(oneside,
           "Plots show regions of (p1,p2) where LNCP is below (yellow), above (red) and near (orange) ",
           "Plots show regions of (p1,p2) where CP is above (yellow), below (red), and near (orange) "),
    ifelse(oneside,alpha/2,1-alpha),
    "\n",
    sep=""),
    line = 5*res.factor,
    outer = TRUE,
    cex=res.factor*textsize)
  }
  mtext(side=2,
        outer=TRUE,
        text = "DNCP\n for individual PSPs",
        cex=textsize*0.8*res.factor,
        at = 1/(2*rows),
        line=0.5*res.factor)
  mtext(side=2,
        outer=TRUE,
        text = "Location index\n for individual PSPs",
        cex=textsize*0.8*res.factor,
        at = 1/(2*rows) + 1/rows,
#        at = (1 - 1/(2*rows) - 2/rows),
        line=0.5*res.factor)
  if (smoothed == TRUE) {
    mtext(side=2,
        outer=TRUE,
        text = paste("Contour plot of \n'Moving Average' ",
                     ifelse(oneside,
                            ifelse(sided=="R","RN","LN"),""),
                     "CP",
                     sep=""),
        cex=textsize*0.8*res.factor,
        at=(1 - 1/(2*rows) - 1/rows),
        line=0.5*res.factor)
  }
  mtext(side=2,
        outer=TRUE,
        text=paste("Plot of ",ifelse(oneside,
                                     ifelse(sided=="R","RN","LN"),""),
                   "CP\n for individual PSPs",
                   sep=""),
        cex=textsize*0.8*res.factor,
        at=(1 - 1/(2*rows)),
        line=0.5*res.factor)
  mtext(side = 3,
        outer = TRUE,
        text = latex2exp::TeX(paste0("\\textbf{",contrast,
                                     "\nN = ", nums,
                                     "\n\u03b1 = ", format(alpha, nsmall=2),
                                     "\n\u03D5 = ", format(par3, nsmall=2),
                                     "}")),
        cex = textsize*0.8*res.factor,
        at = -0.01,
        line = 2*res.factor
        )

  dev.off()
}



