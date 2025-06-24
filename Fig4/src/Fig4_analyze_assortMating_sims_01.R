library(data.table)
library(ggplot2)
library(cowplot)

setwd("~/Dropbox/PENN_STATE_PROJECTS/Matt_admixture_review/assort_mating_slim_code/")

sim = function(rho = 0.99, pop_size = 2000, admix_prop = 0.5) {
  system(paste0("slim -d rho=", rho, " -d total_pop_size=", pop_size," -d admix_prop=", admix_prop, " assortative_mating_ancestry_proportion.slim"))
  ancestr = fread(input = "./output/sim_ancestry_sample_out.csv")
  return(ancestr)
}

rho = 0.9
admix_prop = 0.5
highAssortMating = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 2000, admix_prop = admix_prop))

rho = 0
admix_prop = 0.5
noAssortMating = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 2000, admix_prop = admix_prop))

rho = 0.9
admix_prop = 0.1
highAssortMating_lowAdmixProp = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 2000, admix_prop = admix_prop))

rho = 0
admix_prop = 0.1
noAssortMating_lowAdmixProp = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 2000, admix_prop = admix_prop))

d = cbind(N = 2000, rbind(highAssortMating, noAssortMating, highAssortMating_lowAdmixProp, noAssortMating_lowAdmixProp))


# Also test larger pop size (N = 10000, i.e. less drift)

rho = 0.9
admix_prop = 0.5
highAssortMating = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 10000, admix_prop = admix_prop))

rho = 0
admix_prop = 0.5
noAssortMating = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 10000, admix_prop = admix_prop))

rho = 0.9
admix_prop = 0.1
highAssortMating_lowAdmixProp = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 10000, admix_prop = admix_prop))

rho = 0
admix_prop = 0.1
noAssortMating_lowAdmixProp = cbind(rho = rho, admix_prop = admix_prop, sim(rho, pop_size = 10000, admix_prop = admix_prop))

d = rbind(d, cbind(N = 10000, rbind(highAssortMating, noAssortMating, highAssortMating_lowAdmixProp, noAssortMating_lowAdmixProp)))

# Plots

d[, assortMating := ifelse(rho == 0.0, "Random mating", "High assortative mating")]
d[, admixtureProportion := ifelse(admix_prop == 0.5, "Admixture proportion = 0.5", "Admixture proportion = 0.1")]

d[, assortMating := factor(assortMating, levels = c("Random mating", "High assortative mating"))]
d[, admixtureProportion := factor(admixtureProportion, levels = c("Admixture proportion = 0.5", "Admixture proportion = 0.1"))]


#ggplot(d) + geom_point(aes(Generation, Ancestry_Proportion), alpha = 0.1) + facet_grid(admixtureProportion ~ assortMating) + theme_bw() 


# Plot for paper (N=2000):

pdf(file = "assort_mating_and_admixture.pdf", width = 5, height = 5)
  ggplot(d[Generation < 51]) + geom_point(aes(Generation, Ancestry_Proportion), alpha = 0.1, size = 0.8) + facet_grid(admixtureProportion ~ assortMating) + theme_bw() + 
    labs(x = "Generations after admixture", y = "Individual ancestry proportion") + 
    geom_hline(data = data.frame(admixtureProportion = c("Admixture proportion = 0.5", "Admixture proportion = 0.1"), value = c(0.5, 0.1)), aes(yintercept = value), col="red", lty = 2)
dev.off()


pdf(file = "assort_mating_and_admixture_wideFormat.pdf", width = 9, height = 2.5)
  ggplot(d[Generation < 51]) + geom_point(aes(Generation, Ancestry_Proportion), alpha = 0.1, size = 0.8) + facet_grid(. ~ admixtureProportion + assortMating) + theme_bw() + 
    labs(x = "Generations after admixture", y = "Individual ancestry proportions") + 
    geom_hline(data = data.frame(admixtureProportion = c("Admixture proportion = 0.5", "Admixture proportion = 0.1"), value = c(0.5, 0.1)), aes(yintercept = value), col="red", lty = 2) +
    theme(strip.text.x = element_blank())
dev.off()




# Same plots, but contrasting N=2000 with N=10000 (maybe SI?):

pdf(file = "assort_mating_and_admixture_contrastPopSize.pdf", width = 5, height = 5)
ggplot(d[Generation < 51]) + geom_point(aes(Generation, Ancestry_Proportion, col = factor(N)), alpha = 0.6, size = 0.8) + facet_grid(admixtureProportion ~ assortMating) + theme_bw() + 
  labs(x = "Generations after admixture", y = "Individual ancestry proportion", col = "Population size") + 
  geom_hline(data = data.frame(admixtureProportion = c("Admixture proportion = 0.5", "Admixture proportion = 0.1"), value = c(0.5, 0.1)), aes(yintercept = value), col="red", lty = 2) +
  theme(legend.position = c(0.3, 0.85), legend.background = element_rect(fill = alpha('white', 0.6)))
dev.off()




# Taking samples with different sample size (n = 1, 5, 20, 50, 100), replicated 20 times each:

sampleSizes = rep(c(1,5,20,50,100), each = 20)

d_averages = d[, .(samplesize = sampleSizes, Average_Ancestry_Proportion = sapply(sampleSizes, function(n) mean(sample(size = n, x = Ancestry_Proportion, replace = T)))) 
  , by = .(N, Generation, admixtureProportion, assortMating)]

pdf(file = "assort_mating_and_admixture_when_taking_samples_and_averaging.pdf", width = 5, height = 5)
ggplot(d_averages[Generation < 51 & N == 10000 & samplesize < 50]) + geom_point(aes(Generation, Average_Ancestry_Proportion, col = factor(samplesize)), alpha = 1, size = 0.8) + facet_grid(admixtureProportion ~ assortMating) + theme_bw() + 
  labs(x = "Generations after admixture", y = "Individual ancestry proportion", col = "Sample size") + 
  geom_hline(data = data.frame(admixtureProportion = c("Admixture proportion = 0.5", "Admixture proportion = 0.1"), value = c(0.5, 0.1)), aes(yintercept = value), col="red", lty = 2) +
  theme(legend.position = c(0.3, 0.85), legend.background = element_rect(fill = alpha('white', 0.6)))
dev.off()

