# Pipeline de Bioinformática para Análise de Redes de Co-Expressão Gênica

Este repositório contém o pipeline de bioinformática desenvolvido para a análise de redes de co-expressão gênica usando dados públicos de RNA-seq, com foco em Spodoptera frugiperda, uma importante praga agrícola.

## 🧬 Descrição do Projeto

Com o crescimento do uso de dados públicos de RNA-seq, este pipeline foi criado para:

Obter e filtrar dados do banco SRA (NCBI).
Analisar expressão gênica com pipelines robustos, como nf-core/rnaseq.
Predizer redes de co-expressão gênica usando WGCNA.
O objetivo é explorar os genes associados ao desenvolvimento de Spodoptera frugiperda, contribuindo para pesquisas de controle de pragas.

## 📋 Metodologia

O pipeline é composto pelas seguintes etapas:

1. Download dos dados do SRA e conversão para formato FASTQ:
Obtenção dos dados públicos de RNA-seq e conversão para o formato FASTQ. (NCBI_fastq.sh)

2. Otimização de processamento através de divisão em grupos:
Os dados são organizados em grupos para otimizar a execução e o uso de recursos. (Split_run_size.sh)

3. Execução em ciclo do pipeline nf-core/rnaseq:
Processamento dos dados em ciclos, incluindo alinhamento, quantificação da expressão gênica e controle de qualidade. (NfcoreRnaseqExec.sh)

4. Execução da análise de WGCNA:
Análise das redes de co-expressão gênica para identificar padrões de expressão em diferentes tecidos e fases de vida. (WGCNA.R)

 ## 👨‍💻 Autor

David Daniel Ferreira dos Santos

Universidade: PUC Goiás  
Orientadora: Profª Dra. Mariana Pires de Campos Telles  
Co-Orientadora: Profª Dra. Renata de Oliveira Dias (UFG)  
