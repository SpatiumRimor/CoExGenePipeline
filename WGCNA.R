# Descrição:
# Este script realiza uma análise de expressão gênica utilizando dados de RNA-Seq. Ele inclui as seguintes etapas:
# 1. Instalação e carregamento das bibliotecas necessárias para a análise.
# 2. Leitura dos dados de expressão gênica e metadados.
# 3. Filtragem e limpeza dos dados, incluindo a remoção de genes com contagens zero em todas as amostras.
# 4. Normalização dos dados de expressão utilizando o DESeq2.
# 5. Análise de variância dos genes e filtragem com base no quantil de 95%.
# 6. Geração de gráficos, como violin plots para visualizar a distribuição da expressão gênica por tecido.
# 7. Construção de redes de co-expressão gênica usando WGCNA e identificação de módulos de co-expressão.
# 8. Geração de dendrogramas e visualizações dos módulos.



library(tidyverse)
library(magrittr)  
library(WGCNA)
library(DESeq2)
library(genefilter)

# Carregar os dados de expressão gênica e metadados
data <- readr::read_delim("/media/hd9/star_usage/david-sfru/R/AnalisesFinal/TPMTecidos.tsv",  delim = "\t")
metadata <- readr::read_tsv("/media/hd9/star_usage/david-sfru/R/AnalisesFinal/Tecidosids.txt")

# Verificar as primeiras linhas dos metadados e dados de expressão
head(metadata)
data[1:5, 1:10]

# Renomear a primeira coluna para 'GeneId'
names(data)[1] <- "GeneId"

# Obter os nomes das amostras de expressão (removendo a coluna 'GeneId')
expression_samples <- colnames(data)[-1]

# Identificar amostras faltantes nos metadados
missing_samples <- setdiff(expression_samples, metadata$RUN)
print(missing_samples)

# Identificar amostras comuns entre dados de expressão e metadados
common_samples <- intersect(expression_samples, metadata$RUN)

# Filtrar os dados de expressão para manter apenas as amostras comuns
data_filtered <- data %>% select(GeneId, all_of(common_samples))

# Verificar duplicatas na coluna RUN dos metadados
duplicated_runs <- duplicated(metadata$RUN)
if (any(duplicated_runs)) {
  duplicated_entries <- metadata$RUN[duplicated_runs]
  stop(paste("Metadados contêm RUN duplicados:", paste(unique(duplicated_entries), collapse = ", ")))
} else {
  print("Nenhuma duplicata encontrada na coluna RUN.")
}

# Filtrar e ordenar os metadados para corresponder à ordem das amostras nos dados de expressão
metadata_ordered <- metadata %>%
  filter(RUN %in% common_samples) %>%
  arrange(match(RUN, common_samples))

# Verificar se o ordenamento dos metadados está correto
if (!all(metadata_ordered$RUN == common_samples)) {
  stop("Erro ao ordenar metadados. Verifique se há duplicatas ou inconsistências.")
} else {
  print("Metadados ordenados corretamente.")
}

# Criar o meta_df final, alinhado com os dados de expressão filtrados
meta_df <- data.frame(
  Sample = metadata_ordered$RUN,      # As amostras
  Tissue = metadata_ordered$Tecido    # Os respectivos tecidos
)

# Verificar as primeiras linhas do meta_df
head(meta_df)

# Transformar os dados de expressão para formato longo e mesclar com os metadados
mdata <- data_filtered %>%
  pivot_longer(cols = -GeneId, names_to = "Sample", values_to = "Expression") %>%
  left_join(meta_df, by = "Sample") %>%
  mutate(Tissue = as.factor(Tissue))  # Garantir que os tecidos sejam tratados como fatores

# Ajustar o tamanho da imagem e a resolução
png("plot_output_by_tissue_high_quality.png", width = 35, height = 25, units = "in", res = 600)

# Plotar os dados de expressão com agrupamento por tecido
p <- mdata %>%
  ggplot(aes(x = Sample, y = Expression)) +
  geom_violin() +                                   # Gráfico de violino para mostrar a distribuição dos dados de expressão
  geom_point(alpha = 0.2) +                          # Adicionar pontos individuais para cada amostra, com transparência
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 0.5),  # Ajustar o alinhamento do texto do eixo X
    axis.text.y = element_text(size = 15),                                      # Ajustar o tamanho do texto do eixo Y
    axis.title.x = element_text(size = 30, margin = margin(t = 15)),            # Aumentar o distanciamento da legenda do eixo X
    axis.title.y = element_text(size = 30, margin = margin(r = 15)),            # Aumentar o distanciamento da legenda do eixo Y
    strip.text = element_text(size = 25),                                       # Aumentar o tamanho dos títulos dos grupos
    panel.spacing = unit(1, "lines")                                            # Aumentar o espaçamento entre os gráficos
  ) +
  labs(x = "Samples", y = "RNA Seq Counts") +
  facet_wrap(~ Tissue, scales = "free_x", ncol = 2)  # Facetar o gráfico por tecido, com escalas independentes para cada eixo X

# Salvar o gráfico
print(p)
dev.off()

# Converter os dados de expressão filtrados em uma matriz para análises posteriores
de_input_filtered <- as.matrix(data_filtered[, -1])
rownames(de_input_filtered) <- data_filtered$GeneId

# Verificar quantos genes têm todas as contagens iguais a zero
all_zero_genes <- rowSums(de_input_filtered == 0) == ncol(de_input_filtered)

# Remover genes com todas as contagens iguais a zero
de_input_filtered <- de_input_filtered[!all_zero_genes, ]

# Verificar a dimensão após o filtro
dim(de_input_filtered)

# Adicionar um pseudocount de 1 para evitar problemas com contagens iguais a zero
de_input_filtered_pseudocount <- de_input_filtered + 1

# Criar o DESeqDataSet com a matriz filtrada
dds <- DESeqDataSetFromMatrix(
  countData = round(de_input_filtered_pseudocount),
  colData = meta_df,
  design = ~ Tissue
)

# Executar o DESeq para realizar a normalização dos dados
dds <- DESeq(dds)

# Aplicar a transformação de variância estabilizada para os dados normalizados
vsd <- varianceStabilizingTransformation(dds)

# Extrair os dados transformados para análise posterior
wpn_vsd <- assay(vsd)

# Verificar a dimensão dos dados transformados
dim(wpn_vsd)

# Calcular a variância por gene nos dados transformados
rv_wpn <- rowVars(wpn_vsd)

# Definir o quantil de 95% para filtrar os genes com alta variância
q95_wpn <- quantile(rv_wpn, 0.95)

# Manter apenas os genes cuja variância é maior que o quantil de 95%
expr_normalized <- wpn_vsd[rv_wpn > q95_wpn, ]

# Verificar a dimensão após o filtro
dim(expr_normalized)

# Definir o nome do arquivo PNG, tamanho da imagem e resolução
png("violin_plot_expression.png", width = 10, height = 6, units = "in", res = 300)

# Mesclar os dados de expressão normalizada com os metadados
expr_normalized_df <- data.frame(expr_normalized) %>%
  mutate(Gene_id = row.names(expr_normalized)) %>%
  pivot_longer(-Gene_id, names_to = "Sample", values_to = "Expression") %>%
  left_join(meta_df, by = c("Sample" = "Sample"))  # Mesclar com os metadados

# Gerar o gráfico violin plot usando o tecido no eixo X
png("violin_plot_by_tissue.png", width = 10, height = 6, units = "in", res = 300)

expr_normalized_df %>%
  ggplot(aes(x = Tissue, y = Expression)) +  # Alterar para que o eixo X mostre os tecidos
  geom_violin() +
  geom_point(alpha = 0.2) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90),  # Girar os nomes dos tecidos para facilitar a leitura
    plot.title = element_text(face = "bold", margin = margin(b = 10)),  # Título principal em negrito e afastado
    axis.title.x = element_text(face = "bold", margin = margin(t = 15)),  # Título do eixo X em negrito e afastado
    axis.title.y = element_text(face = "bold", margin = margin(r = 10))   # Título do eixo Y em negrito e afastado
  ) +
  ylim(0, NA) +  # Definir o limite inferior do gráfico como zero
  labs(
    title = "Normalized and 95% Quantile Expression by Tissue",
    x = "Tissue",
    y = "Normalized Expression"
  )

dev.off()  # Salvar o gráfico como PNG

# Transpor a matriz de expressão para que as linhas sejam amostras e as colunas sejam genes
input_mat <- t(expr_normalized)

# Verificar a estrutura da matriz transposta
dim(input_mat)  # Deve retornar [número de amostras] x [número de genes filtrados]

# Permitir multi-threading (opcional) para acelerar a análise
allowWGCNAThreads(12)

# Escolher uma faixa de valores de soft-threshold para testar
powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))

# Analisar a topologia da rede para escolher o soft threshold
sft <- pickSoftThreshold(
  input_mat,
  powerVector = powers,
  verbose = 5
)

# Plotar os resultados da escolha do soft threshold
png("soft_threshold_plots.png", width = 10, height = 6, units = "in", res = 300)
# Definir layout para dois gráficos lado a lado
par(mfrow = c(1, 2))
cex1 <- 0.9  # Tamanho do texto

# Plotar o gráfico de Scale Free Topology Fit
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.56, col = "red")  # Linha de referência para o valor R² = 0.90

# Plotar o gráfico de Mean Connectivity
plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  labels = powers,
  cex = cex1,
  col = "red"
)

dev.off()

# Definir o soft-threshold escolhido (exemplo: power 9)
picked_power <- 9

# Evitar conflitos de namespace com a função 'cor'
temp_cor <- cor       # Salvar a função original 'cor'
cor <- WGCNA::cor     # Usar a função 'cor' do pacote WGCNA

# Construir a rede de co-expressão e identificar módulos
netwk <- blockwiseModules(
  input_mat,
  power = picked_power,
  networkType = "signed",
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minModuleSize = 30,
  maxBlockSize = 4000,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  saveTOMs = TRUE,
  saveTOMFileBase = "ER",
  numericLabels = TRUE,
  verbose = 3
)

# Restaurar a função 'cor' original
cor <- temp_cor

# Converter as etiquetas dos módulos em cores para o gráfico
mergedColors <- labels2colors(netwk$colors)

# Definir o nome do arquivo PNG, tamanho da imagem e resolução
png("dendrogram_modules.png", width = 15, height = 10, units = "in", res = 300)

# Gerar o dendrograma com as cores dos módulos
plotDendroAndColors(
  netwk$dendrograms[[1]],  # Primeiro dendrograma dos blocos
  mergedColors[netwk$blockGenes[[1]]],  # Cores correspondentes aos módulos dos genes
  "Module colors",  # Título do gráfico
  dendroLabels = FALSE,  # Não mostrar labels no dendrograma
  hang = 0.03,  # Controlar o alinhamento das folhas do dendrograma
  addGuide = TRUE,  # Adicionar uma linha-guia para facilitar a visualização
  guideHang = 0.05  # Controlar o comprimento da linha-guia
)

# Finalizar o gráfico e salvar o arquivo
dev.off()

# Gerar o dataframe com genes e seus respectivos módulos
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Salvar a tabela em um arquivo
write_delim(module_df, file = "gene_modules.txt", delim = "\t")

# Obter os Eigengenes por módulo
MEs0 <- moduleEigengenes(input_mat, mergedColors)$eigengenes

# Reordenar os módulos para agrupar módulos similares
MEs0 <- orderMEs(MEs0)
module_order <- names(MEs0) %>% gsub("ME", "", .)

# Incluir os tratamentos e tecidos no dataframe de eigengenes
MEs0$treatment <- row.names(MEs0)
MEs0$Tecido <- metadata_ordered$Tecido  # Adicionar coluna de tecidos usando 'Tecido'

# Transformar os dados de eigengenes para o formato longo
mME <- MEs0 %>%
  pivot_longer(-c(treatment, Tecido)) %>%
  mutate(
    name = gsub("ME", "", name),  # Remover o "ME" dos nomes dos módulos
    name = factor(name, levels = module_order)  # Definir a ordem dos módulos
  )

# Criar o gráfico com os nomes dos tecidos no eixo X
plot <- mME %>%
  ggplot(aes(x = Tecido, y = name, fill = value)) +  # Usar Tecido no eixo X
  geom_tile() +
  theme_bw() +
  scale_fill_gradient2(
    low = "blue",
    high = "red",
    mid = "white",
    midpoint = 0,
    limit = c(-1, 1)
  ) +
  theme(
    axis.text.x = element_text(angle = 90, size = 8),  # Girar e ajustar o tamanho do texto no eixo X
    axis.text.y = element_text(size = 10)  # Ajustar o tamanho do texto no eixo Y
  ) +
  labs(title = "Module-Trait Relationships", y = "Modules", x = "Tissues", fill = "corr")

# Salvar o gráfico como PNG
ggsave("module_trait_relationships_by_tissue.png", plot = plot, width = 10, height = 8, dpi = 300)

# Gerar o dataframe com genes e seus respectivos módulos
module_df <- data.frame(
  gene_id = names(netwk$colors),
  colors = labels2colors(netwk$colors)
)

# Filtrar os dados normalizados para incluir apenas os genes presentes em module_df
submod <- module_df  # Aqui você pode filtrar para módulos específicos, se necessário
subexpr <- expr_normalized[submod$gene_id, ]  # Filtrar os genes normalizados

# Gerar o dataframe com as expressões normalizadas para os genes nos módulos
submod_df <- data.frame(subexpr) %>%
  mutate(gene_id = row.names(.)) %>%
  pivot_longer(-gene_id, names_to = "Sample", values_to = "Expression") %>%
  mutate(module = submod$colors[match(gene_id, submod$gene_id)])

# Verificar se a coluna 'Tissue' existe em 'submod_df' e, se não, adicionar a partir dos metadados
if (!"Tissue" %in% colnames(submod_df)) {
  submod_df <- submod_df %>%
    left_join(meta_df, by = "Sample")  # Supondo que meta_df contém as colunas 'Sample' e 'Tissue'
}

# Adicionar coluna 'LifeStage' com base no nome do tecido
submod_df <- submod_df %>%
  mutate(
    LifeStage = case_when(
      grepl("Larva", Tissue) ~ "Larva",
      grepl("Pupa", Tissue) ~ "Pupa",
      grepl("Adult", Tissue) ~ "Adult",
      TRUE ~ "Other"
    )
  )

# Ajustar o mapeamento de cores dos módulos para tons mais escuros e corretos
module_color_map <- c(
  "blue" = "#1F78B4",        # Azul mais escuro
  "brown" = "#8B4513",       # Marrom escuro
  "green" = "#228B22",       # Verde escuro
  "grey" = "#696969",        # Cinza escuro
  "red" = "#B22222",         # Vermelho escuro
  "turquoise" = "#00CED1",   # Turquesa vibrante
  "yellow" = "#FFD700"       # Amarelo mais escuro e intenso
)

# Definir a ordem dos estágios de vida e ordenar o fator 'Tissue'
submod_df <- submod_df %>%
  mutate(
    LifeStage = factor(LifeStage, levels = c("Larva", "Pupa", "Adult")),
    Tissue = factor(Tissue, levels = unique(Tissue[order(LifeStage)]))
  )

# Gerar o gráfico de expressão normalizada por módulo e tecido, ordenado por estágio de vida
png("expression_by_tissue_life_stage_ordered.png", width = 12, height = 8, units = "in", res = 300)

submod_df %>%
  ggplot(aes(x = Tissue, y = Expression, group = gene_id)) +
  geom_line(aes(color = module), alpha = 0.2) +  # Linhas representando a expressão normalizada de cada gene
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, margin = margin(t = 5)),
    strip.text = element_text(size = 10)
  ) +
  scale_x_discrete(expand = c(0, 0)) +  # Remover espaço extra no eixo X
  scale_color_manual(values = module_color_map) +  # Aplicar o mapeamento de cores personalizado
  facet_grid(rows = vars(module)) +  # Facetar o gráfico por módulo
  labs(
    x = "Tissue (Ordered by Life Stage)",
    y = "Normalized Expression",
    title = "Normalized Expression by Tissue and Module (Ordered by Life Stage)"
  )

dev.off()
# Gráfico para mediada expressão normalizada 

# Calcular a média da expressão normalizada por tecido e módulo
avg_expression_df <- submod_df %>%
  group_by(Tissue, module) %>%
  summarize(MeanExpression = mean(Expression, na.rm = TRUE), .groups = "drop")

# Criar o gráfico png("mean_expression_with_bottom_space.png", width = 14, height = 10, units = "in", res = 300)

avg_expression_df %>%
  ggplot(aes(x = Tissue, y = MeanExpression, group = module)) +
  geom_line(aes(color = module), size = 0.7) +  # Linhas representando a média da expressão normalizada
  geom_point(aes(color = module), size = 2) +  # Pontos para destacar valores médios
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, margin = margin(t = 5)), # Ajustar rotação e tamanho do texto
    axis.title.x = element_text(margin = margin(t = 10)), # Espaçamento para o eixo X
    axis.title.y = element_text(margin = margin(r = 10)), # Espaçamento para o eixo Y
    strip.text = element_text(size = 12, face = "bold"), # Títulos das facetas em negrito e maiores
    plot.margin = margin(t = 30, r = 10, b = 20, l = 10), # Ajustar borda inferior
    panel.spacing = unit(1.2, "lines") # Aumentar o espaçamento entre facetas
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Espaçamento extra no eixo X
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.1))) +  # Espaçamento extra no eixo Y (inferior e superior)
  scale_color_manual(values = module_color_map) +  # Aplicar o mapeamento de cores personalizado
  facet_grid(rows = vars(module)) +  # Facetar o gráfico por módulo
  labs(
    x = "Tissue (Ordered by Life Stage)",
    y = "Mean Normalized Expression",
    title = "Mean Normalized Expression by Tissue and Module (Ordered by Life Stage)" 
  )

dev.off()
