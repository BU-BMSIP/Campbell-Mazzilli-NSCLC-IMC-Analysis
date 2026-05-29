args<-commandArgs(TRUE)
# Pixel Scaling, 99.95%ile method

#Set up directory:
fileDir <- arg[1] 
allFiles <- list.files(fileDir, pattern = ".*ROI.*")

outputDir <- arg[2]


for(ix in 1:length(allFiles)){
    imageTxt <- allFiles[ix]

    txtFile <- read.table(paste0(fileDir, imageTxt), sep = "\t", header = TRUE, check.names = FALSE)

    totalIntensity <- rowSums(txtFile[,7:47])
    totalIntensityThreshold <- quantile(totalIntensity, 0.999)

    #Step 1 - Turn all overly "hot" pixels to 0
    txtFile[which(totalIntensity > totalIntensityThreshold), c(7:47)] <- 0

    #Step 2 - Check for pixels with high intensity on mutually exclusive markers 
    # (Some artifacts have high co-expression of markers not expected to be co-expressing)
    markersToCheck <- c("CD56", "EpCAM", "CD68", "CD117", "CD138", "SMA")
    combnMarkersToCheck <- combn(markersToCheck, 2)

    rmIx <- c()
    for(x in 1:ncol(combnMarkersToCheck)){
        marker1 <- paste0(combnMarkersToCheck[1,x], "\\b")
        marker2 <- paste0(combnMarkersToCheck[2,x], "\\b")
        markerIx1 <- grep(marker1, colnames(txtFile))
        markerIx2 <- grep(marker2, colnames(txtFile))
        intersectionIx <- intersect(which(txtFile[,markerIx1] > quantile(txtFile[,markerIx1], 0.999)),
                                    which(txtFile[,markerIx2] > quantile(txtFile[,markerIx2], 0.999)))
        rmIx <- c(rmIx, intersectionIx)
    }

    # Keep membrane and DNA as is
    txtFile[unique(rmIx), c(7:41,44,45,47)] <- 0
    # 0-100 scaling of intensity
    intensityValuesAdjust <- apply(txtFile[,c(7:41,44,45,47)], 2, function(x){
        return(100 * (x-median(x))/(quantile(x, 0.9995)-median(x)))
        # return(100 * (x-min(x))/(quantile(x, 0.99999)-min(x))) #Used for separate run
    })
    intensityValuesAdjust <- apply(intensityValuesAdjust, 2, function(x) pmax(x, 0))
    intensityValuesAdjust <- apply(intensityValuesAdjust, 2, function(x) pmin(x, 100))
    txtFile[,c(7:41,44, 45, 47)] <- intensityValuesAdjust
    
    
    newtx <- paste0(outputDir, imageTxt)
    #Write into new txt file
    write.table(txtFile, file = newtx, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)
}
