gc()
options (future.globals.maxSize = 2200 * 1024^4) # 设置环境大小限制，稍大于提示信息要求的就行
plan("multisession", workers = 10 ) 
plan()



col =                 c("B"                                     = "#1F78B4",
                        "CD4T"                                  = "#c3d34c",
                        "CD8T"                                  = "#7cd49d",
                        "cDC1"                                  = "#ccc296",
                        "cDC2"                                  = "#d85091",
                        "Cholangiocyte"                         = "#7d5a3b",
                        "Fibroblast"                            = "#38503a",
                        "Gallbladder sinusoid endothelial"      = "#A6CEE3",
                        "Hepatocyte"                            = "#d098ab",
                        "Macrophage"                            = "#3d2130",
                        "Magakaryocyte"                         = "#46de65",
                        "Mast cell"                             = "#cc6e48",
                        "Monocyte"                              = "#ff316c",
                        "Mucosal invariant T"                   = "#746900",
                        "Neutrophil"                            = "#9e71bf",
                        "NK"                                    = "#8f3e7f",
                        "pDC"                                   = "#587f3a",
                        "Plasma cell"                           = "#3f2b71",
                        # "LPM:Extraembryonic mesoderm"             = "#d58d40",
                        # "LPM:Foregut mesenchyme"                  = "#983b40",
                        # "LPM:Gut mesenchyme"                      = "#79bec6",
                        # "LPM:Hepatic mesenchyme"                  = "#663fc6",
                        # "LPM:Proepicardium"                       = "#5b6180",
                        # "LPM:Renal pericytes and mesangial cells" = "#d94e38",
                        # "LPM:Somatic mesoderm"                    = "#7d9cd9",
                        # "LPM:Splanchnic mesoderm"                 = "#ca47cc",
                        # "Neuromesodermal progenitors"             = "#cb566d",
                        # "LPM:Gonad progenitor cells"              = "#42222e",
                        # "LPM:Lung mesenchyme"                     = "#6dd251",
                        # "LPM:Meninges"                            = "#6ad39e",
                        # "LPM:Vascular smooth muscle cells"        = "#c1bf5a",
                        "unknow"                   = "#3a5440")

col              = c("B"                       = "#ca47cc",
                    "CD4T"                     = "#277fb8",
                    "CD8T"                     = "#549da3",
                    "cDC1"                     = "#c6b598",
                    "cDC2"                     = "#ee8e46",
                    "Cholangiocyte"            = "#c1bf5a",
                    "Fibroblast"               = "#9f7ada",
                    "Gallbladder sinusoid endothelial"= "#6ad39e",
                    "Hepatocyte"               = "#dc4989",
                    "Macrophage"               = "#f9a341",
                    "Magakaryocyte"            = "#8f3e7f",
                    "Mast cell"                = "#d58d40",
                    "Monocyte"                 = "#b05a28",
                    "Mucosal invariant T"      = "#a4cde1",
                    "Neutrophil"               = "#79bec6",
                    "NK"                       = "#67a4cc",
                    "pDC"                      = "#7f5d35",
                    "Plasma cell"              = "#42222e",
                    "unknow"                   = "#d399af"
                    # "Granulosa cells"                       = "#516aca",
                    # "Proepicardium"                         = "#5b6180",
                    # "Mesothelial cells"                     = "#3a5440"
                    )
col <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f"
         ,"#8bc96d","#4dae47","#5c9e43","#b79973","#f38989"
         ,"#ec5051","#e32427","#ef6a45","#f9b769","#f9a341"
         ,"#f48521","#ee8e46","#d4a6a8","#af93c4","#8660a8"
         ,"#815e99","#c6b598","#f6f28f","#d4a55b","#b05a28")

col <- c("#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47"
               ,"#5c9e43","#b79973","#f38989","#ec5051","#e32427","#ef6a45","#f9b769"
               ,"#f9a341","#f48521","#ee8e46","#d4a6a8","#af93c4","#8660a8","#815e99"
               ,"#c6b598","#f6f28f","#d4a55b","#b05a28")
