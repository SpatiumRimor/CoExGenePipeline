#!/bin/bash

# Script: NCBI_fastq.sh
# Descrição: Este script realiza o download de arquivos SRA e os converte para o formato FASTQ, verificando se os arquivos já foram baixados e atualizando as listas correspondentes.
# Uso: ./NCBI_fastq.sh


# Verifica se os arquivos 'baixados.txt' e 'SraAccList_novos.txt' existem
if [[ ! -f baixados.txt || ! -f SraAccList_novos.txt ]]; then
    echo "Criando arquivos iniciais 'baixados.txt' e 'SraAccList_novos.txt'."

    # Cria o arquivo 'baixados.txt' com a lista de diretórios existentes
    for dir in */ ; do
        echo "${dir%/}"
    done | sort | uniq > baixados.txt

    # Cria o arquivo 'SraAccList_novos.txt' com os novos itens que ainda não foram baixados
    grep -vxFf baixados.txt SraAccList.txt > SraAccList_novos.txt
else
    echo "Arquivos 'baixados.txt' e 'SraAccList_novos.txt' já existem. Pulando criação inicial."
fi

# Função para verificar se um arquivo está vazio
check_empty_file() {
    if [[ ! -s $1 ]]; then
        echo "Nenhum arquivo para processar. Todos os arquivos foram baixados com sucesso."
        exit 0
    fi
}

# Função para processar cada linha do arquivo (cada ID SRA)
process_line() {
    local sra_id=$1
    local retry=0

    while true; do
        echo
        echo "Baixando $sra_id"
        echo

        # Baixando arquivo com prefetch e capturando toda a saída
        echo "Iniciando prefetch para $sra_id..."
        prefetch $sra_id --max-size 1000000000000000000000000000000000000
        local prefetch_exit_code=$?

        # Verifica se o comando prefetch foi bem-sucedido
        if [[ $prefetch_exit_code -ne 0 ]]; then
            # Se houver erro de serviço externo, tenta novamente
            if prefetch $sra_id 2>&1 | grep -q "Failed to call external services."; then
                echo "Falha ao chamar serviços externos durante prefetch. Tentativa: $((retry+1)) para $sra_id."
                echo
                sleep 10
                ((retry++))
                continue
            else
                # Se houver outro erro inesperado, retorna falha
                echo "Erro inesperado durante prefetch"
                echo
                return 1
            fi
        fi

        # Convertendo arquivo SRA para fastq com fasterq-dump
        fasterq-dump $sra_id
        if [[ $? -eq 3 ]]; then
            # Se o arquivo já existir, realiza comandos adicionais e atualiza a lista
            echo "Erro de arquivo já existente. Executando comandos adicionais."
            echo # Adiciona uma linha em branco

            # Atualiza a lista de arquivos baixados
            update_downloaded_list

            return 0
        fi

        # Se tudo deu certo, sai do loop
        break
    done
}

# Função para atualizar a lista de arquivos baixados
update_downloaded_list() {
    # Atualiza 'baixados.txt' com a lista de diretórios baixados
    for dir in */ ; do
        echo "${dir%/}"
    done | sort | uniq > baixados.txt

    # Atualiza 'SraAccList_novos.txt' com os IDs que ainda precisam ser baixados
    grep -vxFf baixados.txt SraAccList.txt > SraAccList_novos.txt
}

# Loop externo que continuará até que todos os arquivos sejam baixados
while true; do
    # Verifica se há arquivos a serem processados
    check_empty_file "SraAccList_novos.txt"

    # Loop para ler o arquivo de entrada 'SraAccList_novos.txt' linha por linha
    while IFS=$'\t' read -r linha1 linha2 linha3 linha4 linha5; do
        need_retry=false
        while true; do
            # Processa cada linha, que corresponde a um ID SRA
            process_line "$linha1"
            
            # Se o process_line retornar 1, significa que houve falha e precisa tentar novamente
            if [[ $? -eq 1 ]]; then
                echo "Reiniciando tentativa para o arquivo $linha1."
                echo
                need_retry=true
                sleep 10
                continue
            fi

            # Se o process_line retornar 0, significa que o download foi bem-sucedido
            break
        done
    done < SraAccList_novos.txt

    # Atualiza a lista de arquivos baixados após tentar todos os IDs
    update_downloaded_list

    # Verifica novamente se há arquivos a serem processados
    if [[ ! -s SraAccList_novos.txt ]]; then
        echo "Todos os arquivos foram baixados com sucesso."
        break
    fi

    # Se houver falhas, tenta baixar novamente
    if $need_retry; then
        echo "Tentando novamente baixar arquivos que falharam..."
        echo
    fi
done
