#!/bin/bash
set -e  # Faz o script parar em caso de erro

# Autor: David Daniel
# Descrição: Este script automatiza o processamento e empacotamento de runs de dados de RNA-Seq usando o pipeline nf-core/rnaseq. 
# Ele realiza as seguintes etapas principais:
# 1. Verifica se cada run já foi processada e registrada em um log.
# 2. Exclui arquivos temporários ('work') para liberar espaço.
# 3. Cria um pacote tar.gz dos dados processados para armazenamento eficiente.
# 4. Processa os dados de RNA-Seq com o Nextflow, gerando relatórios de qualidade com o MultiQC.
# 5. Registra as runs processadas em um log.

# Definição de Cores para Melhorar a Estética da Saída
GREEN='[0;32m'  # Cor verde para mensagens de sucesso
YELLOW='[1;33m'  # Cor amarela para avisos
RED='[0;31m'     # Cor vermelha para erros
BLUE='[0;34m'    # Cor azul para informações
NC='[0m'         # Sem cor

# Funções para Exibir Mensagens com Cores
function echo_info() {
    echo -e "${BLUE}$1${NC}"
}

function echo_success() {
    echo -e "${GREEN}$1${NC}"
}

function echo_warning() {
    echo -e "${YELLOW}$1${NC}"
}

function echo_error() {
    echo -e "${RED}$1${NC}"
}

# Função para Verificar Integridade do tar.gz usando tar -tzf
function verify_tar_integrity() {
    local tar_file="$1"
    
    # Verifica se o arquivo não está vazio
    if [ ! -s "$tar_file" ]; then
        echo_error "Arquivo $tar_file está vazio."
        return 1
    fi
    
    # Tenta listar o conteúdo do tar.gz
    if ! tar -tzf "$tar_file" > /dev/null 2>&1; then
        echo_error "Arquivo $tar_file está corrompido ou inválido."
        return 1
    fi
    
    echo_success "Arquivo $tar_file está íntegro."
    return 0
}

# Definições de Diretórios e Arquivos
BASE_DIR=$(pwd)  # Define o diretório base como o diretório atual
SAMPLE_SHEET_ORIGINAL="/media/hd9/star_usage/david-sfru/Requisitos/sample_sheet.csv"  # Caminho do arquivo de sample sheet original
LOG_FILE="$BASE_DIR/processed_runs.log"  # Arquivo de log para registrar runs processadas
DEST_DIR="/media/lgbio-nas1/davidsantos/TCC/Star-Nfcore/Analises"  # Diretório de destino para os arquivos empacotados

# Verificar se o diretório base existe
if [ ! -d "$BASE_DIR" ]; then
    echo_error "Erro: O diretório base $BASE_DIR não existe."
    exit 1
fi

# Verificar se o diretório de destino existe
if [ ! -d "$DEST_DIR" ]; then
    echo_error "Erro: O diretório de destino $DEST_DIR não existe."
    exit 1
fi

# Criar o arquivo de log se não existir
touch "$LOG_FILE"

# Navegar para o diretório base
cd "$BASE_DIR" || { echo_error "Erro: Não foi possível acessar $BASE_DIR."; exit 1; }

# Gerar a lista de Runs
Runs=$(ls Run*.txt 2>/dev/null | grep -oP '\d+' | sort -n)  # Lista todos os arquivos Run*.txt, extrai os números e os ordena

# Verificar se existem Runs para processar
if [ -z "$Runs" ]; then
    echo_error "Erro: Nenhum arquivo Run*.txt encontrado em $BASE_DIR."
    exit 1
fi

# Função para Verificar se uma Run já foi Processada
function is_run_processed() {
    local run_number="$1"
    grep -qw "Run$run_number" "$LOG_FILE"  # Verifica se a run já está registrada no log
}

# Loop para cada Run
for Run in $Runs; do
    echo "----------------------------------------"
    echo_info "Processando Run$Run"
    echo "----------------------------------------"

    # Verificar se a Run já foi processada
    if is_run_processed "$Run"; then
        echo_success "Run$Run já foi processada anteriormente. Pulando para a próxima Run..."
        continue
    fi

    # Definir variáveis
    RunDir="$BASE_DIR/Run$Run"  # Diretório da run atual
    OutputDir="$RunDir/star-rsemN$Run"  # Diretório de saída para a run atual
    MultiqcReport="$OutputDir/multiqc/star_rsem/multiqc_report.html"  # Caminho do relatório MultiQC
    RunTxt="$BASE_DIR/Run$Run.txt"  # Caminho do arquivo de entrada da run
    TarFileDest="$DEST_DIR/Run$Run.tar.gz"  # Caminho final do arquivo tar.gz no destino
    TempTarFileDest="$BASE_DIR/Run$Run.tar.gz.tmp"  # Caminho temporário do arquivo tar.gz no diretório base

    # Verificar se o arquivo Run$Run.txt existe
    if [ ! -f "$RunTxt" ]; then
        echo_error "Erro: O arquivo $RunTxt não foi encontrado."
        continue
    fi

    # Verificar se o arquivo tar.gz já existe e está íntegro no destino
    if [ -f "$TarFileDest" ]; then
        echo_info "Verificando a integridade do arquivo tar existente no destino: $TarFileDest"
        if verify_tar_integrity "$TarFileDest"; then
            echo_success "Run$Run já foi empacotada e está íntegra. Adicionando ao log e pulando para a próxima Run..."
            echo "Run$Run" >> "$LOG_FILE"
            continue
        else
            echo_warning "Run$Run possui um arquivo tar existente no destino, mas ele está corrompido ou incompleto. Recriando o tar.gz..."
            # Remove o tar.gz corrompido no destino
            rm -f "$TarFileDest"
            if [ $? -ne 0 ]; then
                echo_error "Erro ao remover o arquivo tar corrompido: $TarFileDest"
                continue
            fi
        fi
    fi

    # Verificar se o relatório MultiQC já existe
    if [ -f "$MultiqcReport" ]; then
        echo_success "Relatório MultiQC para Run$Run já existe em $MultiqcReport. Empacotando a Run..."

        # ============================================
        # EXCLUI A PASTA 'work' para otimização de espaço
        # ============================================
        WorkDir="$RunDir/work"
        if [ -d "$WorkDir" ]; then
            echo_info "Pasta 'work' encontrada em $RunDir. Excluindo..."
            rm -rf "$WorkDir"
            if [ $? -eq 0 ]; then
                echo_success "Pasta 'work' excluída com sucesso."
            else
                echo_error "Erro ao excluir a pasta 'work' em $RunDir."
                continue  # Pular para a próxima Run se não for possível excluir
            fi
        else
            echo_info "Pasta 'work' não encontrada em $RunDir."
        fi
        # ============================================

        # Empacotar a Run
        echo_info "Criando arquivo tar temporário em $BASE_DIR para Run$Run..."
        tar -czvf "$TempTarFileDest" -C "$BASE_DIR" "Run$Run"  # Cria um arquivo tar.gz temporário contendo a pasta da run
        if [ $? -ne 0 ]; then
            echo_error "Erro: Falha ao criar o arquivo tar temporário: $TempTarFileDest"
            # Remove o arquivo temporário se a criação falhar
            rm -f "$TempTarFileDest"
            continue
        fi

        echo_success "Arquivo tar temporário criado com sucesso: $TempTarFileDest"

        # Verificar a integridade do arquivo tar temporário usando tar -tzf
        echo_info "Verificando a integridade do arquivo tar temporário: $TempTarFileDest"
        if verify_tar_integrity "$TempTarFileDest"; then
            echo_success "Arquivo tar temporário $TempTarFileDest está válido."
        else
            echo_error "Erro: O arquivo tar temporário $TempTarFileDest está corrompido."
            # Remove o arquivo temporário
            rm -f "$TempTarFileDest"
            continue
        fi

        # Mover o arquivo tar temporário para o destino
        echo_info "Movendo o arquivo tar temporário para o diretório final: $TarFileDest"
        mv "$TempTarFileDest" "$TarFileDest"  # Move o arquivo tar.gz temporário para o diretório de destino
        if [ $? -ne 0 ]; then
            echo_error "Erro ao mover o arquivo tar para o diretório final: $TarFileDest"
            # Remove o arquivo temporário caso o movimento falhe
            rm -f "$TempTarFileDest"
            continue
        fi

        echo_success "Arquivo tar movido com sucesso para $TarFileDest."

        # Excluir a pasta Run após a compressão bem-sucedida
        echo_info "Excluindo a pasta Run$Run..."
        rm -rf "$RunDir"  # Exclui o diretório da run após a compressão bem-sucedida
        if [ $? -eq 0 ]; then
            echo_success "Pasta Run$Run excluída com sucesso."
        else
            echo_error "Erro ao excluir a pasta Run$Run."
        fi

        # Registrar a Run como processada no log (opcional)
        echo "Run$Run" >> "$LOG_FILE"  # Adiciona a run ao log de processados
        echo_success "Run$Run marcada como concluída no log."

        continue
    fi

    # Se o relatório MultiQC não existe, processar com Nextflow
    echo_info "Relatório MultiQC para Run$Run não encontrado. Iniciando processamento com Nextflow."

    # Criar diretório para o Run (se não existir)
    mkdir -p "$RunDir"

    # Mudar para o diretório do Run
    cd "$RunDir" || { echo_error "Erro: Não foi possível acessar $RunDir."; continue; }

    # Executar o comando grep e adicionar o cabeçalho
    echo_info "Gerando sample_sheet.csv para Run$Run..."
    grep -F -f "$RunTxt" "$SAMPLE_SHEET_ORIGINAL" | \
    sed "1i\sample,fastq_1,fastq_2,strandedness" > sample_sheet.csv  # Cria o arquivo sample_sheet.csv com o cabeçalho adequado

    # Verificar se sample_sheet.csv foi criado corretamente
    if [ ! -s "sample_sheet.csv" ]; then
        echo_error "Erro: sample_sheet.csv está vazio ou não foi criado corretamente para Run$Run."
        cd "$BASE_DIR"
        continue
    fi

    # Criar links simbólicos para os arquivos listados em Run$Run.txt
    echo_info "Criando links simbólicos para arquivos listados em $RunTxt..."
    while IFS= read -r line; do
        SOURCE_FILE="/media/lgbio-nas1/davidsantos/TCC/RNASeq/$line"  # Caminho do arquivo de origem
        if [ -f "$SOURCE_FILE" ]; then
            ln -sf "$SOURCE_FILE" .  # Cria um link simbólico para cada arquivo listado
        else
            echo_warning "Aviso: O arquivo $SOURCE_FILE não existe. Link simbólico não criado."
        fi
    done < "$RunTxt"

    # Criar links simbólicos para os arquivos necessários
    echo_info "Criando links simbólicos para arquivos de requisitos..."
    ln -sf /media/hd9/star_usage/david-sfru/Requisitos/sfru/ .
    ln -sf /media/hd9/star_usage/david-sfru/Requisitos/Sfru.gtf.gz .
    ln -sf /media/hd9/star_usage/david-sfru/Requisitos/Sfru.fna .
    ln -sf /media/hd9/star_usage/david-sfru/Requisitos/rsem .

    # Executar o comando Nextflow com Feedback Visual
    echo_info "Iniciando execução do Nextflow para Run$Run..."
    nextflow run nf-core/rnaseq \
        -c /media/hd12-Renata/lgbio_guest/Nf-Core/Time.config \ # Arquivo de configuração para aumentar tempo de execução para etapas longas 
        --input sample_sheet.csv \
        --aligner star_rsem \
        --star_index sfru \
        --rsem_index rsem \
        --fasta Sfru.fna \
        --gtf Sfru.gtf.gz \
        --outdir "$OutputDir" \
        --skip_biotype_qc \
        --remove_ribo_rna \
        --max_cpus 20 \
        -profile singularity \
        --max_time 200h \
        -resume  # Resume o pipeline em caso de erro

    # Após a execução, verificar se o relatório MultiQC foi gerado
    if [ -f "$MultiqcReport" ]; then
        echo_success "Processamento de Run$Run finalizado com sucesso. Relatório MultiQC encontrado em $MultiqcReport."
    else
        echo_error "Erro: O relatório MultiQC para Run$Run não foi encontrado em $MultiqcReport."
        cd "$BASE_DIR"
        continue
    fi

    # Excluir a pasta 'work' da Run atual antes de compactar
    echo_info "Excluindo a pasta 'work' de Run$Run..."
    WORK_DIR="$RunDir/work"

    if [ -d "$WORK_DIR" ]; then
        rm -rf "$WORK_DIR"  # Exclui a pasta 'work' para liberar espaço
        if [ $? -eq 0 ]; then
            echo_success "Pasta 'work' excluída com sucesso para Run$Run."
        else
            echo_error "Erro: Falha ao excluir a pasta 'work' para Run$Run."
            cd "$BASE_DIR"
            continue
        fi
    else
        echo_warning "Pasta 'work' não encontrada em $RunDir. Nenhuma exclusão necessária."
    fi

    # Voltar para o diretório base antes de empacotar
    cd "$BASE_DIR" || { echo_error "Erro: Não foi possível retornar para $BASE_DIR."; exit 1; }

    # Empacotamento do Run
    echo_info "Iniciando o empacotamento do Run$Run..."

    # Criar arquivo tar.gz temporário no BASE_DIR
    echo_info "Criando arquivo tar temporário em $BASE_DIR para Run$Run..."
    tar -czvf "$TempTarFileDest" -C "$BASE_DIR" "Run$Run"  # Cria um arquivo tar.gz contendo a pasta da run
    if [ $? -ne 0 ]; then
        echo_error "Erro: Falha ao criar o arquivo tar temporário: $TempTarFileDest"
        # Remove o arquivo temporário se a criação falhar
        rm -f "$TempTarFileDest"
        continue
    fi

    echo_success "Arquivo tar temporário criado com sucesso: $TempTarFileDest"

    # Verificar a integridade do arquivo tar temporário usando tar -tzf
    echo_info "Verificando a integridade do arquivo tar temporário: $TempTarFileDest"
    if verify_tar_integrity "$TempTarFileDest"; then
        echo_success "Arquivo tar temporário $TempTarFileDest está válido."
    else
        echo_error "Erro: O arquivo tar temporário $TempTarFileDest está corrompido."
        # Remove o arquivo temporário
        rm -f "$TempTarFileDest"
        continue
    fi

    # Mover o arquivo tar temporário para o destino
    echo_info "Movendo o arquivo tar temporário para o diretório final: $TarFileDest"
    mv "$TempTarFileDest" "$TarFileDest"  # Move o arquivo tar.gz para o diretório de destino
    if [ $? -ne 0 ];then
        echo_error "Erro ao mover o arquivo tar para o diretório final: $TarFileDest"
        # Remove o arquivo temporário caso o movimento falhe
        rm -f "$TempTarFileDest"
        continue
    fi

    echo_success "Arquivo tar movido com sucesso para $TarFileDest."

    # Excluir a pasta Run após a compressão bem-sucedida
    echo_info "Excluindo a pasta Run$Run..."
    rm -rf "$RunDir"  # Exclui o diretório da run após a compactação
    if [ $? -eq 0 ]; then
        echo_success "Pasta Run$Run excluída com sucesso."
    else
        echo_error "Erro ao excluir a pasta Run$Run."
    fi

    # Registrar a Run como processada no log (opcional)
    echo "Run$Run" >> "$LOG_FILE"  # Adiciona a run ao log de processados
    echo_success "Run$Run marcada como concluída no log."

done

echo "----------------------------------------"
echo_success "Todos os Runs foram processados com sucesso."
echo "----------------------------------------"
