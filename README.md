Pipeline de Bioinformática para Análise de Redes de Co-Expressão Gênica

Este repositório contém o pipeline de bioinformática desenvolvido para a análise de redes de co-expressão gênica usando dados públicos de RNA-seq, com foco em Spodoptera frugiperda, uma importante praga agrícola.

🧬 Descrição do Projeto
Com o crescimento do uso de dados públicos de RNA-seq, este pipeline foi criado para:

Obter e filtrar dados do banco SRA (NCBI).
Analisar expressão gênica com pipelines robustos, como nf-core/rnaseq.
Predizer redes de co-expressão gênica usando WGCNA.
O objetivo é explorar os genes associados ao desenvolvimento de Spodoptera frugiperda, contribuindo para pesquisas de controle de pragas.

📋 Metodologia
O pipeline é composto pelas seguintes etapas:

1. Obtenção dos Dados
Busca no NCBI-SRA por dados de RNA-seq relacionados a S. frugiperda.
Automação do download via script Bash usando prefetch e fasterq-dump.
2. Pré-Processamento dos Dados
Remoção de leituras de baixa qualidade com CutAdapt e Trimmomatic.
Filtragem de sequências contaminantes com critérios baseados no conteúdo de GC e cobertura genômica.
3. Análise de Expressão Gênica
Utilização do pipeline nf-core/rnaseq para alinhamento e quantificação dos dados de expressão gênica.
Ferramentas utilizadas:
Alinhamento: STAR
Quantificação: RSEM
Controle de Qualidade: FastQC, SortMeRNA, MultiQC
4. Análise de Redes de Co-Expressão
Predição de módulos de co-expressão gênica com WGCNA.
Identificação de genes centrais e módulos associados a tecidos e fases de desenvolvimento.

👨‍💻 Autor
David Daniel Ferreira dos Santos

Universidade: PUC Goiás
Orientadora: Profª Dra. Mariana Pires de Campos Telles
Co-Orientadora: Profª Dra. Renata de Oliveira Dias (UFG)
