# Cargamos pacman
library( "pacman" )

# cargamos mas paquetes
p_load( "vroom",
        "purrr",
        "tidyr",
        "dplyr",
        "viridis",
        "pheatmap",
        "ggplot2" )

# Cargamos todos los tsv en un solo dataframe
filenames <- list.files( pattern="*.tsv" )

alineamientos <- map_df( filenames,
                         vroom ) %>% 
  separate( data = .,
            col = referenceID,
            into = c( "reference",
                      "protein_ref" ),
            sep = "_" ) %>% 
  separate( data = .,
            col = query_seqID,
            into = c( "query_species",
                      "protein_query" ),
            sep = "_" ) %>% 
  filter( reference != query_species )

# Hacemos heatmap
# creamos el heatmap
elplot <- ggplot( data = alineamientos,
                  mapping = aes( x = protein_query,
                                 y = query_species,
                                 fill = identity ) ) +
  geom_tile( color = "white",
             size = 2 ) +
  labs( title = "Global alignments",
        subtitle = paste( "vs the reference:",
                          alineamientos$reference %>% unique( )  ) ) +
  # scale_fill_gradient( low = "tomato",
  # high = "gray" )
  scale_x_discrete( position = "top" ) +
  scale_fill_viridis_c( option = "mako",
                        limits = c( 50, 100) ) +
  theme_bw( ) +
  theme( legend.title = element_blank( ),
         plot.title = element_text( hjust = 0.5 ),
         plot.subtitle = element_text( hjust = 0.5 ) )

#vis
elplot

# guardamosplot
ggsave( filename = "heatmap_normal.png" ,
        width = 7,
        height = 7,
        dpi = 600 )

# Pasamos a formato wide
ancho <- pivot_wider( data = alineamientos,
                      id_cols = c(reference, query_species),
                      names_from = protein_query,
                      values_from = identity  )

# intentamos hacer clustering por pheatmap
# nos quedamos solo con la matriz de numeros
matriz <- ancho %>% 
  select( 3:ncol( . ) )

# le ponemos nombres a las filas de la matriz
row.names( matriz ) <- ancho$query_species

# abrimos el device para guardar el pheatmap
pdf( file = "pheatmap_clustering.pdf",
     width = 7, height = 7 )

# hacemos pheatmap :)
pheatmap( mat = matriz,
          border_color = "gray30",
          cellwidth = 50,
          scale = "none",
          main = paste( "Global alignments vs", alineamientos$reference %>% unique( ) ),
          color = mako( 50 ),
          breaks = seq( from = 50,
                        to = 100,
                        by = 1 ),
          cutree_rows = 2 )

# cerramos el pdf
dev.off( )

# guardamos la tabla ancha
write.csv( x = ancho,
           file = "identity.csv",
           quote = FALSE,
           row.names = FALSE )
