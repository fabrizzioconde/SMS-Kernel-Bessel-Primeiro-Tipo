# Função que compara as duas listas e retorna os SNPs causais selecionados e não selecionados
selecionar_snp_causais <- function(selected_snps, causal_snps) {
  # Encontrar os SNPs causais que foram selecionados (interseção entre as listas)
  snps_causais_selecionados <- intersect(selected_snps, causal_snps)
  
  # Encontrar os SNPs causais que NÃO foram selecionados (elementos de causal_snps ausentes em selected_snps)
  snps_causais_nao_selecionados <- setdiff(causal_snps, selected_snps)
  
  # Retorna uma lista com os dois resultados
  return(list(
    "SNPs_causais_selecionados" = snps_causais_selecionados,
    "SNPs_causais_nao_selecionados" = snps_causais_nao_selecionados
  ))
}

for (i in 1:6)
{
  # Exemplo de uso:
  # Lista de SNPs selecionados pelo algoritmo SMS
  selected_snps <- snps_selec_corte[[i]]
  
  # Lista de SNPs causais
  causal_snps <- c("SNP1", "SNP10","SNP20", "SNP30", "SNP40", "SNP50", "SNP60")
  
  # Chamada da função
  resultado <- selecionar_snp_causais(selected_snps, causal_snps)
  
  # Impressão dos resultados
  cat("SNPs causais selecionados:", resultado$SNPs_causais_selecionados, "\n")
  cat("SNPs causais não selecionados:", resultado$SNPs_causais_nao_selecionados, "\n")
}
