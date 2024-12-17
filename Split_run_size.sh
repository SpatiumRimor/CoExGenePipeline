#!/bin/bash

# Script: split_run_size.sh
# Descrição: Divide uma lista de arquivos fastq.gz em múltiplos arquivos de execução (Run) com base no tamanho total especificado.
# Uso: ./split_run_size.sh fastqSizesRuns.txt XXXGB X

# Verifica se os argumentos necessários foram fornecidos
if [ "$#" -ne 3 ]; then
    echo "Uso: $0 fastqSizesRuns.txt TAMANHO_MAXIMO NUMERO_INICIAL (exemplo: 200GB 5)"
    exit 1
fi

INPUT_FILE="$1"
MAX_SIZE="$2"
START_RUN_NUMBER="$3"

# Verifica se o arquivo de entrada existe
echo "Verificando se o arquivo de entrada existe: $INPUT_FILE"
if [ ! -f "$INPUT_FILE" ]; then
    echo "Erro: Arquivo de entrada $INPUT_FILE não encontrado."
    exit 1
fi

# Função para converter tamanho para bytes
# Esta função recebe uma string representando um tamanho (e.g., "200GB") e converte para bytes
size_to_bytes() {
    local size_str="$1"
    size_str=$(echo "$size_str" | tr -d "'\"")  # Remove aspas
    size_str=$(echo "$size_str" | sed 's/,/./g')  # Substitui vírgulas por pontos
    number=$(echo "$size_str" | sed -E 's/^([0-9.]+)[KMGT]?B?$/\1/')  # Extrai a parte numérica
    unit=$(echo "$size_str" | sed -E 's/^[0-9.]+([KMGT]?B?)$/\1/')  # Extrai a unidade
    case "$unit" in
        K|KB|k|kb) bytes=$(awk "BEGIN {print $number * 1000}") ;;
        M|MB|m|mb) bytes=$(awk "BEGIN {print $number * 1000 * 1000}") ;;
        G|GB|g|gb) bytes=$(awk "BEGIN {print $number * 1000 * 1000 * 1000}") ;;
        T|TB|t|tb) bytes=$(awk "BEGIN {print $number * 1000 * 1000 * 1000 * 1000}") ;;
        B|b) bytes=$(awk "BEGIN {print $number}") ;;
        *) bytes=0 ;;
    esac
    echo "$bytes"
}

# Função para converter bytes para formato legível
# Esta função converte um valor em bytes para um formato mais legível (e.g., "200GB")
bytes_to_human() {
    local bytes=$1
    local unit="B"
    local size=$bytes
    if (( bytes >= 1000**4 )); then
        size=$(awk "BEGIN {printf \"%.2f\", $bytes/1000000000000}")
        unit="TB"
    elif (( bytes >= 1000**3 )); then
        size=$(awk "BEGIN {printf \"%.2f\", $bytes/1000000000}")
        unit="GB"
    elif (( bytes >= 1000**2 )); then
        size=$(awk "BEGIN {printf \"%.2f\", $bytes/1000000}")
        unit="MB"
    elif (( bytes >= 1000 )); then
        size=$(awk "BEGIN {printf \"%.2f\", $bytes/1000}")
        unit="KB"
    fi
    echo "${size}${unit}"
}

# Converte o TAMANHO_MAXIMO para bytes
echo "Convertendo TAMANHO_MAXIMO para bytes: $MAX_SIZE"
MAX_SIZE_BYTES=$(size_to_bytes "$MAX_SIZE")
if [ "$MAX_SIZE_BYTES" -eq 0 ]; then
    echo "Erro: Formato de TAMANHO_MAXIMO inválido. Use formatos como 200GB, 100MB, etc."
    exit 1
fi

# Inicializa variáveis
echo "Inicializando variáveis"
run_number=$START_RUN_NUMBER  # Começa do número especificado
current_size=0
output_file="Run${run_number}.txt"
> "$output_file"  # Cria ou esvazia o arquivo de saída inicial

declare -A run_sizes  # Declaração de um array associativo para armazenar os tamanhos dos arquivos de execução

prev_id=""
group_files=()
group_size=0

# Função para adicionar um grupo de arquivos ao run atual
# Verifica se o tamanho do grupo ultrapassa o limite e cria um novo run, se necessário
add_group_to_run() {
    local group_size=$1
    local -n group_files_arr=$2  # Referência para o array de arquivos

    # Verifica se o grupo atual excede o tamanho máximo permitido
    if (( current_size + group_size > MAX_SIZE_BYTES )); then
        if (( current_size == 0 )); then
            # Caso em que o tamanho do grupo é maior que o permitido e não há um run atual
            echo "O tamanho do grupo é maior que o tamanho máximo permitido e nenhum run atual existente."
        else
            # Finaliza o run atual e cria um novo
            run_sizes["Run${run_number}.txt"]=$current_size
            size_human=$(bytes_to_human "$current_size")
            echo "Run${run_number}.txt criado com tamanho: $size_human"
            run_number=$((run_number + 1))
            output_file="Run${run_number}.txt"
            > "$output_file"
            current_size=0
        fi
    fi

    # Adiciona cada arquivo do grupo ao run atual
    for f in "${group_files_arr[@]}"; do
        echo "Adicionando arquivo ao run atual: $f"
        echo "$f" >> "$output_file"
    done
    current_size=$((current_size + group_size))
}

line_count=0

echo "Iniciando processamento do arquivo de entrada"
# Loop para processar cada linha do arquivo de entrada
while IFS= read -r line || [ -n "$line" ]; do
    line_count=$((line_count + 1))
    echo "Processando linha $line_count: $line"
    if [ -z "$line" ]; then
        echo "Linha vazia encontrada, ignorando"
        continue
    fi

    # Extrai o nome do arquivo e o tamanho da linha atual
    filename=$(echo "$line" | awk '{print $1}')
    size_str=$(echo "$line" | awk '{print $2}' | tr -d "'\"" | sed 's/,/./g')

    # Converte o tamanho para bytes
    size_bytes=$(size_to_bytes "$size_str")

    if [ "$size_bytes" -eq 0 ]; then
        echo "Aviso: Formato de tamanho '$size_str' para o arquivo '$filename' inválido. Ignorando."
        continue
    fi

    # Extrai o ID do arquivo, removendo sufixos como "_1" ou "_2" e a extensão .fastq.gz
    id=$(echo "$filename" | sed -E 's/(_[12])?\.fastq\.gz$//')

    echo "ID do arquivo: $id"
    # Verifica se o ID é igual ao do arquivo anterior, indicando que faz parte do mesmo grupo
    if [ "$id" == "$prev_id" ]; then
        echo "ID igual ao anterior ($prev_id), adicionando ao grupo atual"
        group_files+=("$filename")
        group_size=$((group_size + size_bytes))
        add_group_to_run "$group_size" group_files
        group_files=()
        group_size=0
        prev_id=""
    else
        # Caso o ID seja diferente, finaliza o grupo anterior (se houver) e inicia um novo
        if [ "${#group_files[@]}" -gt 0 ]; then
            echo "Novo ID encontrado, adicionando grupo atual ao run e iniciando novo grupo"
            add_group_to_run "$group_size" group_files
            group_files=()
            group_size=0
        fi
        group_files=("$filename")
        group_size="$size_bytes"
        prev_id="$id"
    fi
done < "$INPUT_FILE"

# Adiciona o último grupo ao run, se houver
if [ "${#group_files[@]}" -gt 0 ]; then
    echo "Adicionando grupo final ao run"
    add_group_to_run "$group_size" group_files
fi

# Finaliza o último run
run_sizes["Run${run_number}.txt"]=$current_size
size_human=$(bytes_to_human "$current_size")
echo "Run${run_number}.txt criado com tamanho: $size_human"
echo "Total de linhas processadas: $line_count"
echo "Divisão concluída em Run${START_RUN_NUMBER}.txt até Run${run_number}.txt"
