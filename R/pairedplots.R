

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

if(FALSE){
plotdata=myarrays
alpha=0.05
nums="30,30"
xlim=c(0,1)
ylim=xlim
methlab="MN"
avg=F
lside=F
lines=F
lines1=F
CIlen=T
crude=F
res.factor=6
colour=T
texsize=1
}

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
	# OScols to be updated for alpha=0.1 (and 0.01?)
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
#		              col2,col3,
                  col2,col2,col3,col3,
		              rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(8),each=1))),
		"0.05" = rev(c(rep("#000000",380),
		               rep(colorRampPalette(c("BLACK",col1))(9)[-1],each=1),
		               col2,col2,col3,col3,
		               rep(colorRampPalette(c(col4,"WHITE"),space="Lab")(8),each=1))),
		"0.01" = rev(c(rep("#000000",1980),
		               rep(colorRampPalette(c("BLACK",col1))(9)[-1],each=1),
#		               col2,col3,
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
#	par(bg="lightgray") # Set colour for missing values
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
#	if (any(is.na(ccol))) {
#	  mmat <- ifelse(is.na(ccol), 1, NA)
#	  image(mmat, axes = FALSE, xlab = "", ylab = "", col = "lightgray", useRaster=TRUE, add = TRUE)
#	}
#	if (any(is.na(rawcpdata))) {
#	  mmat <- ifelse((is.na(rawcpdata) | is.na(ccol)), 1, NA)
#	  image(mmat, axes = FALSE, xlab = "", ylab = "", col = "lightgray", useRaster=TRUE, add = TRUE)
#	}
	if (avg | CIlen) {
# if (FALSE) {
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
#		if (lside) 	abline(coef=c(0.1,1),lwd=0.5*textsize*res.factor)
		if (lside && lines %in% c("both","1")) {
		  lines(x=c(0.5,0.85), y=c(0.6,0.95), lwd = 0.5*textsize*res.factor)
#		}	else if (lside==FALSE & lines %in% c("both","2")) {
		}	else if (lside==FALSE & lines %in% c("RR")) {
		  #		  abline(coef = c(1.4,-1), lwd = 0.5*textsize*res.factor)
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

plotpanel <- function(plotdata,
                      alpha,
                      par3,
#                      nums,
                      sel,
                      oneside = F,
                      plotlab,
                      linesx = F,
                      smoothed = F,
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
#  longlab[longlab=="Tang-cc5"] <- "Tang-cc (\U0263=0.5)"
#  longlab[longlab=="Tang-cc25"] <- "Tang-cc (\U0263=0.25)"
#  longlab[longlab=="Tang-cc125"] <- "Tang-cc (\U0263=0.125)"
  longlab[longlab=="SCAS"] <- "SCASu"
  longlab[longlab=="SCAS-bc"] <- "SCAS"
  longlab[longlab=="MOVER-NJcc5"] <- "MOVER-c (\U0263=0.5)"
  longlab[longlab=="MOVER-NJcc125"] <- "MOVER-c (\U0263=0.125)"
  longlab[longlab=="SCAS-cc5"] <- "SCAS-c (\U0263=0.5)"
  longlab[longlab=="SCAS-cc25"] <- "SCAS-c (\U0263=0.25)"
  longlab[longlab=="SCAS-cc125"] <- "SCAS-c (\U0263=0.125)"
  if (FALSE) {
    longlab[longlab=="SC"] <- "SCAS"
    #	longlab[longlab=="SC"] <- "GNbc"
    longlab[longlab=="SCcc"] <- "SCAS-cc (\U0263=0.5)"
    longlab[longlab=="SCcomp"] <- "SCAS-cc (\U0263=0.25)"
    longlab[longlab=="Wald"] <- "AN"
    longlab[longlab=="AMA"] <- "AMB-approx"
    longlab[longlab=="AMB"] <- "AMB-exact"
    longlab[longlab=="MOVER-E"] <- "MOVER-cc (\U0263=0.5)"
    longlab[longlab=="MOVER-comp"] <- "MOVER-cc (\U0263=0.25)"
    longlab[longlab=="Wa"] <- "AN"
    longlab[longlab=="Woolf"] <- "AN"
    longlab[longlab=="Katz"] <- "AN"
    #	longlab=c("Miettinen-Nurminen","Wallenstein (bc)","Agresti-Caffo","Newcombe-Wilson","Mee","Wallenstein",
    #	"MN (N/N-2)","Wallenstein (bc2)","Wald","Profile Likelihood","Newcombe-Wilson (CC)",
    #	"Wald (CC)","Agresti-Min","Chan-Zhang","Edgeworth Expansion","Jeffreys-Perks","Brown-Li")
  }
  names(longlab) <- names(methods) <- methods
  n.grid <- dim(plotdata[["mastercp"]])[1]
  if (colour == F) collab=".bw" else collab=""

  rows <- 4

  if(fmt=="tiff")  {
    tiff(file = paste(outpath,"_",fmt,"/","summary",
                      ifelse(sided=="L","L",""),plotlab,100*(1-alpha),"_",
                      nums,ifelse(oneside,"os",""),collab,".tiff",
                      sep=""
    ),
    width=600*res.factor,
    height=330*res.factor,
    type="quartz"
    ) #,res=360,compression="none")  ## to create a tiff file ##   changed from 600x360
  } else if(fmt=="png") {
    png(file = paste(outpath,"_","png","/","summary",ifelse(sided=="L","L",""),
                     plotlab,nums,"_",par3,"_",100*(1-alpha),
                     ifelse(oneside,"os",""),collab,".png",
                     sep=""
    ),
    width = (120 * nmeth + 60) * res.factor,
#    height = 12 * 40 * res.factor,
    height = 4*rows * 38 * res.factor,
    type = "quartz"
    ) #,res=360,compression="none")  ## to create a png file ##   changed from 600x360
  } else if(fmt=="eps") {
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
  layout(matrix(c(1:(rows * nmeth), rows * nmeth + c(1, 1:(rows-1))), rows, (nmeth + 1),
                byrow = FALSE),
         widths = c(rep(2, nmeth), 1), heights = rep(4, rows))
#  widths=c(2,2,2,2,1))
  par(oma = res.factor*c(0,3,ifelse(fmt=="xxx",10,7),0), pty='s')
  par(cex.main = res.factor*0.8*textsize, cex.axis=res.factor*0.8*textsize)
  par(mar = res.factor*(c(2,1,1,0.5)+0.1))
#  par(bg = "lightgray") # Set colour for missing values

  #	i <- 1
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
#                        ifelse(oneside,
#                               paste0("\n",ifelse(sided=="R","maxRNCP=","maxLNCP="),
#                                      summaries[i,"maxLNCP"]),
#                               paste0("\n","meanCP=",summaries[i,"meanCP"],
#                                      "\n","minCP=",summaries[i,"minCP"])
#                        ),
            "aveCons=",summaries[i,"pctAveCons"],"%",
            ifelse(any(sel=="SCcc"),
                               paste0(("\n"),"pctCons=",(summaries[i,paste("pctCons",ifelse(oneside,".1side",""),sep="")]),"%"),
                               paste0("\n","ave within ±",
                                      (ifelse(oneside,
                                              format(0.1*alpha/2,scientific=F),
                                              format(0.1*alpha,scientific=F))
                                      ),"=",
                                      (summaries[i,paste("pctAvenear",ifelse(oneside,".1side",""),sep="")]),
                                      "%"
                               )
                        )))
    )

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
  #          "\n","mean DNCP=",
  #          summaries[i,"meanDNCP"]
#            "\n","DNCP bad=",
#            (summaries[i,"pctBad.DNCP"]),
#            "%"
#  "\n","DNCP above ",
#  format(1.2*alpha/2,scientific=F),
#  "=",(summaries[i,"pctBad.DNCP"]),
#  "%"
  #    "\n","1-sided prox=",
        #    (summaries[i,"pctnear.1side"]),
        #    "%" meanDNCP
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
#                          locind = TRUE,
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
            "%"
          )))
  }

  par(pty="m", mar=res.factor*c(3,2,ifelse(fmt=="xxx",2,2),3)+0.1)

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
#           y=(c(180:200)-180)/20,
#           labels=(180:200)/200,
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
           #y=(c(80:100)-80)/20,
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
#      image(y=(1:50-0.5)/50,
#            z=matrix(1951:2000,nrow=1),
#            col=rev(palette)[1951:2000],
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
#      image(y=(1:100-0.5)/100,
#            z=matrix(1901:2000,nrow=1),
#            col=palette[1901:2000],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10, 0.45, 0.55),
#           y=(c(190:200)-190)/10,
#           labels=(190:200)/200,
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
#            col = rev(palettex)[181:200],
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
#      image(y=(1:50-0.5)/50,
#            z=matrix(1951:2000,nrow=1),
#            col=rev(palettex)[1951:2000],
            axes=F,
            ylab="")
      text(x=1.2,
           y = c(0:10/10),
#           labels = c(((5*98):(5*100))/(5*100)),
#           y=(c(195:200)-195)/5,
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
        at=1/(2*rows),
        line=0.5*res.factor)
  mtext(side=2,
        outer=TRUE,
        text = "Location index\n for individual PSPs",
        cex=textsize*0.8*res.factor,
        at = (1 - 1/(2*rows) - 2/rows),
        line=0.5*res.factor)
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
        text = latex2exp::TeX(paste0("\\textbf{",contrast, "\nN=", nums,"}")),
        cex = textsize*0.8*res.factor,
        at = -0.01,
#        adj = 0,
        line = 4.5*res.factor
        )

  dev.off()
}

if (FALSE) {
  # Keep old version
plotpanel_old <- function(plotdata,
                      alpha,
                      psi,
                      nums,
                      sel,
                      oneside,
                      plotlab,
                      linesx = F,
                      CIlen = F,
                      res.factor = 6,
                      fmt = "png",
                      colour = T,
                      textsize = 1,
                      limits = c(0,1),
                      sided = "R"
                      ) {

  #plotdata=myarrays;alpha=0.05;nums="15,15";sel=c("MN","AC","N","WL");oneside=T;plotlab="a"
#  plotdata=myarrays;alpha=0.05; psi=10;nums=1;sel=c("Tang","Tang-bc","MOVER-nw","MOVER-ns");oneside=F;plotlab="a"
  #linesx=T;CIlen=F;res.factor=6;
#  textsize=1
  g <- "gamma"
	summaries <- plotdata[["summaries"]][paste(psi),,paste(100*(1-alpha)),,nums,]
	longlab <- methods <- dimnames(plotdata[["mastercp"]])[[4]]
	if (FALSE) {
	  longlab[longlab=="SC"] <- "SCAS"
	  longlab[longlab=="MOVER-NJc"] <- "MOVER-NJ"
	  #	longlab[longlab=="SC"] <- "GNbc"
	longlab[longlab=="SCcc"] <- "SCAS-cc (\U0263=0.5)"
	longlab[longlab=="SCcomp"] <- "SCAS-cc (\U0263=0.25)"
	longlab[longlab=="Wald"] <- "AN"
	longlab[longlab=="AMA"] <- "AMB-approx"
	longlab[longlab=="AMB"] <- "AMB-exact"
	longlab[longlab=="MOVER-E"] <- "MOVER-cc (\U0263=0.5)"
	longlab[longlab=="MOVER-comp"] <- "MOVER-cc (\U0263=0.25)"
	longlab[longlab=="Wa"] <- "AN"
	longlab[longlab=="Woolf"] <- "AN"
	longlab[longlab=="Katz"] <- "AN"
#	longlab=c("Miettinen-Nurminen","Wallenstein (bc)","Agresti-Caffo","Newcombe-Wilson","Mee","Wallenstein",
#	"MN (N/N-2)","Wallenstein (bc2)","Wald","Profile Likelihood","Newcombe-Wilson (CC)",
#	"Wald (CC)","Agresti-Min","Chan-Zhang","Edgeworth Expansion","Jeffreys-Perks","Brown-Li")
	}
	names(longlab) <- names(methods) <- methods
	n.grid <- dim(plotdata[["mastercp"]])[1]
	if (colour == F) collab=".bw" else collab=""

  if(fmt=="tiff")  {
    tiff(file = paste(outpath,"_",fmt,"/","summary",
                      ifelse(sided=="L","L",""),plotlab,100*(1-alpha),"_",
                    nums,ifelse(oneside,"os",""),collab,".tiff",
                    sep=""
                    ),
         width=600*res.factor,
         height=330*res.factor,
         type="quartz"
         ) #,res=360,compression="none")  ## to create a tiff file ##   changed from 600x360
  } else if(fmt=="png") {
    png(file = paste(outpath,"_","png","/","summary",
                     ifelse(sided=="L","L",""), plotlab, 100*(1-alpha), "_",
                     psi,"_",
                     nums,ifelse(oneside,"os",""),collab,".png",
                     sep=""
                     ),
        width=600*res.factor,
        height=360*res.factor,
        type="quartz"
        ) #,res=360,compression="none")  ## to create a png file ##   changed from 600x360
  } else if(fmt=="eps") {
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

	layout(matrix(c(1,2,3,4,5,6,7,8,9,9), 2, 5,
	              byrow = FALSE),
	       widths=c(2,2,2,2,1))
#alpha.lab=bquote(paste(.(ifelse(oneside,"",paste("\n","aveCons=",summaries[i,"pctAveCons"],"%",sep=""))),.("\n"),.("ave within +/-"),"0.1",alpha))
	par(oma = res.factor*c(0,4,ifelse(fmt=="xxx",10,4),0), pty='s')
	par(cex.main = res.factor*0.8*textsize,cex.axis=res.factor*0.8*textsize)
	par(mar = res.factor*(c(2,1,1,0.5)+0.1))
#	i <- methods["Tang"]
	for(i in methods[sel]) {
	  # Run the unsmoothed plot and capture the output palette
		palette <- CPcontour(plotdata=plotdata,
		                     alpha=alpha,
		                     psi=psi,
		                     nums=nums,
		                     methlab=i,
		                     lside=oneside,
		                     lines=linesx,
		                     CIlen=CIlen,
		                     res.factor=res.factor,
		                     colour=colour,
		                     textsize=textsize,
		                     xlim=limits)
		mtext(side = 3,
		      cex = res.factor*0.8*textsize,
		      line = 1*res.factor,
		      text = (paste(longlab[i],
		              ifelse(oneside,
		                         paste("\n",ifelse(sided=="R","maxRNCP=","maxLNCP="),
		                                       summaries[i,"maxLNCP"],sep=""),
		                         paste("\n","minCP=",summaries[i,"minCP"],sep="")
		                         ),
		              ifelse(any(sel=="SCcc"),
		                     paste0(("\n"),"pctCons=",(summaries[i,paste("pctCons",ifelse(oneside,".1side",""),sep="")]),"%"),
		                     paste0(("\n"),"within ±",(ifelse(oneside,format(0.1*alpha/2,scientific=F),format(0.1*alpha,scientific=F))),"=",(summaries[i,paste("pctnear",ifelse(oneside,".1side",""),sep="")]),"%")
		                     ),
		              sep=""))
		      )

		# Run the moving average smoothed plot
		CPcontour(plotdata = plotdata,
		          alpha = alpha,
		          psi=psi,
		          nums = nums,
		          methlab = i,
		          lside = oneside,
		          avg = T,
		          lines = linesx,
		          CIlen = CIlen,
		          res.factor = res.factor,
		          colour = colour,
		          textsize = textsize,
		          xlim = limits)

		mtext(side=3,
		      cex=res.factor*0.8*textsize,
		      line=0.5*res.factor,
		      text = (paste(
		        ifelse(oneside,
		               paste("\n","meanloc=",summaries[i,"meanlocindex"],sep=""),
		               paste("\n","aveCons=",summaries[i,"pctAveCons"],"%",sep="")
		               ),
		        "\n","ave within ±",
		        ifelse(oneside,
		                format(0.1*alpha/2,scientific=F),
		                format(0.1*alpha,scientific=F)
		               ),
		        "=",(summaries[i,paste("pctAvenear",
		                               ifelse(oneside,".1side",""),
		                               sep="")]),
		        "%",
		        sep=""))
		)
	}


	par(pty="m", mar=res.factor*c(2,2,ifelse(fmt=="xxx",2,1),3)+0.1)

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
			     y=(c(180:200)-180)/20,
			     labels=(180:200)/200,
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
			     y=(c(80:100)-80)/20,
			     labels=(80:100)/100,
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
			image(y=(1:50-0.5)/50,
			      z=matrix(1951:2000,nrow=1),
			      col=rev(palette)[1951:2000],
			      axes=F,
			      ylab="")
			text(x=1.2,
			     y=(c(195:200)-195)/5,
			     labels=(5:0)/200,
			     xpd=T,
			     cex=res.factor*textsize,
			     pos=4)
		} else {
			image(y=(1:100-0.5)/100,
			      z=matrix(1901:2000,nrow=1),
			      col=palette[1901:2000],
			      axes=F,
			      ylab="")
			text(x=1.2,
			     y=(c(190:200)-190)/10,
			     labels=(190:200)/200,
			     xpd=T,
			     cex=res.factor*textsize,
			     pos=4)
		}
	}
	box(lwd=0.6*textsize*res.factor)
	mtext(text=ifelse(oneside,
	                  ifelse(sided=="R","RNCP","LNCP"),"CP"),
	      side=3,
	      at=2,
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
	      text = paste("Contour plot of \n'Moving Average' ",
	                   ifelse(oneside,
	                          ifelse(sided=="R","RN","LN"),""),
	                   "CP",
	                   sep=""),
	      cex=textsize*0.8*res.factor,
	      at=0.25,
	      line=1*res.factor)
	mtext(side=2,
	      outer=TRUE,
	      text=paste("Plot of ",ifelse(oneside,
	                                   ifelse(sided=="R","RN","LN"),""),
	                 "CP\n for individual PSPs",
	                 sep=""),
	      cex=textsize*0.8*res.factor,
	      at=0.75,
	      line=1*res.factor)

	dev.off()
}
}

if (FALSE) {
RRpairteam <- c("Tang", "Tang-bc", "Tang-sc", "Tang-scbc") 	#Paired RR
RR2pairteam <- c("MOVER-w", "MOVER-nw", "MOVER-nj","MOVER-ns") 	#Paired RR
RRccpairteam <- c("Tang-cc5", "Tang-cc125", "Tang-sccc", "MOVER-ccns") 	#Paired RR, cc

RDpairteam <- c("Tango", "Tango-bc", "Tango-sc", "Tango-scbc") 	#Paired RD
RD2pairteam <- c("MOVER-w", "MOVER-nw", "MOVER-nj","MOVER-ns") 	#Paired RD
RDccpairteam <- c("Tango-cc5", "Tango-cc125", "Tango-sccc", "MOVER-ccns") 	#Paired RD, cc
}

if (FALSE) {

for (i in c(0.01, 0.05, 0.1)) {
#  for (i in c(FALSE, TRUE)) {
  for (j in c(2, 10)) {
    plotpanel(plotdata=cparrays_RD20,alpha=i, psi=j,nums="20",limits= plotlim,sel=RDpairteam,oneside=F,plotlab="RDpair",res.factor=res.factor,fmt=fmt,linesx="RD",colour=colour,sided=sided)
  }
}
#}

for (j in c(2, 10)) {
  for (i in c(0.05, 0.1)) {
    #for (j in c(2, 10, 100)) {
    #      for (i in c(0.01, 0.05, 0.1)) {
    #for (i in c(FALSE, TRUE)) {
    plotpanel(plotdata=cparrays_RR40, alpha=i, psi=j, nums="40",
              limits= plotlim,sel=RRpairteam,oneside=F,plotlab="RRpair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
    plotpanel(plotdata=cparrays_RR40, alpha=i, psi=j, nums="40",
              limits= plotlim,sel=RR2pairteam,oneside=F,plotlab="RR2pair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
    plotpanel(plotdata=cparrays_RR40, alpha=i, psi=j, nums="40",
              limits= plotlim,sel=RRccpairteam,oneside=F,plotlab="RRccpair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
  }
}

teamlist <- list(RDpairteam, RD2pairteam, RDccpairteam)
teamlabels <- c("RDpair", "RD2pair", "RDccpair")

teamlist <- list(RRpairteam, RR2pairteam, RRccpairteam)
teamlabels <- c("RRpair", "RR2pair", "RRccpair")

for (j in c(2, 10)) {
  for (i in c(0.1, 0.05)) {
    for (k in 1:3) {
    plotpanel(plotdata=cparrays_RR20, alpha=i, psi=j, nums="20",
              limits= plotlim,sel=teamlist[[k]],oneside=F,plotlab=teamlabels[k],
              res.factor=res.factor,fmt=fmt,linesx=F,colour=colour,sided=sided)
    }
  }
}

dev.off()

for (j in c(2, 10)) {
  for (i in c(0.1)) {
    #for (j in c(2, 10, 100)) {
    #      for (i in c(0.01, 0.05, 0.1)) {
    #for (i in c(FALSE, TRUE)) {
    plotpanel(plotdata=cparrays_RR20, alpha=i, psi=j, nums="20",
              limits= plotlim,sel=RRpairteam,oneside=F,plotlab="RRpair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
    plotpanel(plotdata=cparrays_RR20, alpha=i, psi=j, nums="20",
              limits= plotlim,sel=RR2pairteam,oneside=F,plotlab="RR2pair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
    plotpanel(plotdata=cparrays_RR20, alpha=i, psi=j, nums="20",
              limits= plotlim,sel=RRccpairteam,oneside=F,plotlab="RRccpair",
              res.factor=res.factor,fmt=fmt,linesx="RR",colour=colour,sided=sided)
  }
}


#}

plotpanel(plotdata=cparrays_RR60,alpha=0.05, psi=10,nums="60",limits= plotlim,sel=RRpairteam,oneside=F,plotlab="RRpair",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
dev.off()

plotpanel(plotdata=cparrays_RR40,alpha=0.05, psi=100,nums=1,limits= plotlim,sel=RRpairteam,oneside=F,plotlab="RRpair",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)
plotpanel(plotdata=cparrays_RR40,alpha=0.05, psi=100,nums=1,limits= plotlim,sel=RRpairteam,oneside=T,plotlab="RRpair",res.factor=res.factor,fmt=fmt,linesx=lines,colour=colour,sided=sided)

}

if(FALSE) {


#outpath="/Users/Pete/Main/Contract/ssu2010_003/DiffBin/plots/"
#path="/Users/Pete/Main/Contract/ssu2010_003/DiffBin/R/"
m=15;n=15
m=20;n=10
m=120;n=60
load(file=paste(outpath,"masterarraysx.RR.",m,".",n,".Rdata",sep=""))
p1=as.numeric(dimnames(myarrays[[1]])[[1]])
p2=as.numeric(dimnames(myarrays[[1]])[[2]])
methods=dimnames(myarrays[[1]])[[3]]
nmeth<-length(methods)

nums<-dimnames(myarrays[[2]])[[4]]

nums=paste(m,",",n,sep="")
alpha=0.1
oneside=T
oneside=F
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("N","MN","BLJ","MN/BLJ"),oneside=oneside,plotlab="a")
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("Wall","PL","AC","BLR"),oneside=oneside,plotlab="b")
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("WaldCC","NCC","CZ","AM"),oneside=oneside,plotlab="c")

#nn<-800
#round(diffBinconf.all(0.78*nn,0.27*nn,nn,nn),3)[c("MN","Wbc","BLj","N"),,]


plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("MN","BL","N","Wbc"),oneside=oneside,plotlab="a")

plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("AM","CZ","WaldCC","NCC"),oneside=oneside,plotlab="e")
#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("OR|p2","OR|psi","HW","MN"),oneside=oneside,plotlab="d")
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("Wald","PL","Mee","Wall"),oneside=oneside,plotlab="c")
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("MN","JP","EE","Wbc"),oneside=oneside,plotlab="f")
plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("Mee","Wall","MN2","Wbc2"),oneside=oneside,plotlab="b")

sapply(1:dim(myarrays[[2]])[2],function(i) {
	nums<-dimnames(myarrays[[2]])[[2]][i]
	#plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("WL"),oneside=oneside,plotlab="a")
	plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("MN","Wbc","BL","N"),oneside=oneside,plotlab="a")
	plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("Wald","PL","AC","JP"),oneside=oneside,plotlab="b")
	plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("Mee","Wall","MN2","Wbc2"),oneside=oneside,plotlab="c")
#	plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("AM","CZ","WaldCC","NCC"),oneside=oneside,plotlab="d")
#	plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("OR|p2","OR|psi","HW","MN"),oneside=oneside,plotlab="d")
	}
)

plotpanel(plotdata=myarrays,alpha=alpha,nums=nums,sel=c("MN","AC","N","WL"),oneside=T,plotlab="b")
plotpanel(plotdata=myarrays,alpha=0.05,nums="15,15",sel=c(1:min(4,nmeth)),oneside=F,plotlab=paste("RR","len",sep=""),CIlen=T)

if (nmeth>4) {
	plotpanel(sel=c(5:min(8,nmeth)),oneside=F,plotlab="b")
	plotpanel(sel=c(5:min(8,nmeth)),oneside=T,plotlab="b")
	plotpanel(sel=c(5:min(8,nmeth)),oneside=F,plotlab=paste("b","len",sep=""),CIlen=T)
}


#some stratified plots:
alpha=0.05
oneside=F
nums="15,10"
plotpanel(plotdata=myarrays.strat,alpha=alpha,nums=nums,sel=c("MNCMH","MNinv","MNshrink","Crude"),oneside=oneside,plotlab="sa")
plotpanel(plotdata=myarrays.strat,alpha=alpha,nums=nums,sel=c("BLCMH","WaldShrink","Wald.metabin","Wald.me"),oneside=oneside,plotlab="sb")
plotpanel(plotdata=myarrays.strat,alpha=alpha,nums=nums,sel=c("WLCMH","WLshrink","WLinv","Wald.AC"),oneside=oneside,plotlab="sc")
plotpanel(plotdata=myarrays.strat,alpha=alpha,nums=nums,sel=c("BLCMH","BLshrink","BLinv","Crude"),oneside=oneside,plotlab="sd")

dimnames(myarrays.strat[[2]])

dev.off()

}
