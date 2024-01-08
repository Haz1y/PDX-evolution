
my_gistic2map_noheader = function (X,refBuild = 'hg19') {
    raw_gistic = read_tsv(X)
    raw_gistic$Chromosome = gsub(pattern = 'X', replacement = '23', x = raw_gistic$Chromosome, fixed = TRUE)

    raw_gistic_amp = raw_gistic %>% filter(Type == "Amp") %>% as.data.frame
    raw_gistic_del = raw_gistic %>% filter(Type == "Del") %>% as.data.frame

    ## Change End
    for (i in (2: nrow(raw_gistic_amp)-1)) {
     if ( (raw_gistic_amp$Chromosome[i] == raw_gistic_amp$Chromosome[(i+1)])) { 
        raw_gistic_amp$End[i] = raw_gistic_amp$Start[i+1]
         }
      else {raw_gistic_amp$End[i] = raw_gistic_amp$Start[i]}
        }

    for (i in (2: nrow(raw_gistic_del)-1)) {
         if ( (raw_gistic_del$Chromosome[i] == raw_gistic_del$Chromosome[(i+1)])) { 
            raw_gistic_del$End[i] = raw_gistic_del$Start[i+1]
             }
          else {raw_gistic_del$End[i] = raw_gistic_del$Start[i]}
        }
    
    ## 
    if(refBuild == 'hg19'){
        chrLens = c(249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663,
                    146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540,
                    102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566,
                    155270560)
                    #, 59373566)
    } else if(refBuild == 'hg18'){
        chrLens = c(247249719, 242951149, 199501827, 191273063, 180857866, 170899992,
                    158821424, 146274826, 140273252, 135374737, 134452384, 132349534,
                    114142980, 106368585, 100338915, 88827254, 78774742, 76117153,
                    63811651, 62435964, 46944323, 49691432, 154913754)
                    #, 57772954)
    } else if(refBuild == 'hg38'){ 
        chrLens = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979,
                    159345973, 145138636, 138394717, 133797422, 135086622, 133275309,
                    114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
                    58617616, 64444167, 46709983, 50818468, 156040895)
                    #, 57227415)
    }
    
    chrLabels <- seq_len(23)

    chrTable <- data.table::data.table(chr = as.character(chrLabels),
                               start = as.numeric(chrLens),
                               end = as.numeric(chrLens))
    chrLens <- append(cumsum(chrTable$start),1,after = 0)
    chrTable <- chrTable %>% 
        dplyr::mutate(
            start = chrLens[seq_len(length(chrLens)-1)],
            end = chrLens[2:length(chrLens)]
        )

    chrLabels <- chrTable$chr
    
    ## color 
    chr_bar_color <- lapply(seq_len(length(as.numeric(chrLabels) )), function(i){
        if((i %% 2) == 0){
            return("white")
        }else{
            return("gray90")
        }
    }) %>% unlist()

    chr_text_color <- lapply(seq_len(length(as.numeric(chrLabels))), function(i){
        if((i %% 2) == 0){
            return("gray90")
        }else{
            return("white")
        }
    }) %>% unlist()


    chrTable$color = chr_bar_color
    chr_start_pos <- as.numeric(chrTable$start) 
    names(chr_start_pos) <- paste0("chr",chrTable$chr)

    
    ## rearrange 
    raw_gistic_amp <- raw_gistic_amp %>% 
        dplyr::mutate(
            Update_Start = as.numeric(.data$Start) + chr_start_pos[paste0("chr",.data$Chromosome)],
            Update_End = as.numeric(.data$End) + chr_start_pos[paste0("chr",.data$Chromosome)]
        )

    raw_gistic_del <- raw_gistic_del %>% 
        dplyr::mutate(
            Update_Start = as.numeric(.data$Start) + chr_start_pos[paste0("chr",.data$Chromosome)],
            Update_End = as.numeric(.data$End) + chr_start_pos[paste0("chr",.data$Chromosome)]
        )
    
    p = ggplot() + 
    geom_rect(data = raw_gistic_amp ,
              aes(xmin = Update_Start,
                  xmax = Update_End,
                  ymin = 0, ymax = frequency),
             fill = "#DA3A32",col = "#DA3A32")+
    geom_rect(data = raw_gistic_del ,
              aes(xmin = Update_Start,
                  xmax = Update_End,
                  ymin = -frequency, ymax = 0),
             fill = "#3F529F",col = "#3F529F")+
    geom_segment(aes(x = chrTable$end[seq_len(nrow(chrTable)-1)],
                     xend = chrTable$end[seq_len(nrow(chrTable)-1)],
                     y = -Inf,
                     yend = 1),
                 #linetype = "dotted",
                 color = "black",
                 size = 0.7)+
    geom_segment(aes(x = centromeres.hg19_pos$start[seq_len(nrow(chrTable))],
                 xend = centromeres.hg19_pos$start[seq_len(nrow(chrTable))],
                 y = -Inf,
                 yend = 1),
             linetype = "dotted",
             color = "gray30",
             size = 0.7)+
    theme(
            axis.text.x = element_blank(),
            # axis.text.y.left =  element_text(size = sample.text.size,
            #                                  margin = margin(l = 0),
            #                                  colour = "black"),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(colour = "black", size = 0.5),
            axis.ticks.length.y=unit(0.25, "cm"),
            axis.line.y = element_line(colour = "black", size = 0.5),
            panel.grid =element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            # legend.key.width  = unit(1.1,'lines'),
            # legend.box.margin = margin(l = 0),
            # legend.title = element_text(size = legend.title.size),
            # legend.text = element_text(size = legend.text.size,margin = margin(b = 3))
        )+
        # coord_fixed(ratio = 1) + 
        scale_x_continuous(expand = c(0,0))+scale_y_continuous(limit = c(-0.8,1),position = "right")
        # scale_y_continuous(breaks = Y.text.table$Pos,
        #                    labels = Y.text.table$Tumor_Sample_Barcode)+
        #ggtitle("Copy number alteration profile")+
        theme(plot.title = element_text(size = 13.5, face = "bold", hjust = 0.5, vjust = -2))
   return (p) 
    
}
