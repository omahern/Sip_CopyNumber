sample_data(physeq_rep3_t)

incorp_only_1=incorp_only


glu_asvs=subset(incorp_only_1, id=='batch_glu')$OTU
sub1=subset_taxa(physeq_rep3_t, Strain %in% glu_asvs)
glu_only=subset_samples(sub1, Culture == "Batch" & Treatment=="Single")
glu_melt <- psmelt(glu_only)
glu_melt <- glu_melt %>%
  group_by(Strain, Isotope) %>%                       # .grad = one tube (isotope + replicate)
  mutate(gmax = max(Abundance, na.rm = TRUE),
         rel  = ifelse(gmax > 0, Abundance / gmax, 0)) %>%
  ungroup()


multi_batch_asvs=subset(incorp_only_1, id=='batch_carb')$OTU
sub1=subset_taxa(physeq_rep3_t, Strain %in% multi_batch_asvs)
batch_multi=subset_samples(sub1, Culture == "Batch" & Treatment=="Multi")
batch_multi_melt <- psmelt(batch_multi)
batch_multi_melt <- batch_multi_melt %>%
  group_by(Strain, Isotope) %>%                       # .grad = one tube (isotope + replicate)
  mutate(gmax = max(Abundance, na.rm = TRUE),
         rel  = ifelse(gmax > 0, Abundance / gmax, 0)) %>%
  ungroup()


multi_chemo_asvs=subset(incorp_only_1, id=='chemo_carb')$OTU
sub1=subset_taxa(physeq_rep3_t, Strain %in% multi_chemo_asvs)
chemo_multi=subset_samples(sub1, Culture == "Chemostat" & Treatment=="Multi")
chemo_multi_melt <- psmelt(chemo_multi)

chemo_multi_melt <- chemo_multi_melt %>%
  group_by(Strain, Isotope) %>%                       # .grad = one tube (isotope + replicate)
  mutate(gmax = max(Abundance, na.rm = TRUE),
         rel  = ifelse(gmax > 0, Abundance / gmax, 0)) %>%
  ungroup()


colors_fill=c("C12" = "gray", "C13" = 'firebrick')


library(ggpubr)
ggplot(glu_melt, aes(x=Buoyant_density, y=rel, fill=Isotope, color=Isotope)) +
  geom_point() +
  geom_line() +
  theme_classic2(base_size = 10) +
  facet_wrap(~Strain) +
  scale_fill_manual(values=colors_fill) +
  scale_color_manual(values=colors_fill) +
  ylab("Ratio of Maximum Quantity") +
  xlab("Buoyant Density (g/mL)") +
  theme(legend.position = 'bottom')

library(ggpubr)
ggplot(batch_multi_melt, aes(x=Buoyant_density, y=rel, fill=Isotope, color=Isotope)) +
  geom_point() +
  geom_line() +
  theme_classic2(base_size = 10) +
  facet_wrap(~Strain) +
  scale_fill_manual(values=colors_fill) +
  scale_color_manual(values=colors_fill) +
  ylab("Ratio of Maximum Quantity") +
  xlab("Buoyant Density (g/mL)") +
  theme(legend.position = 'bottom')

library(ggpubr)
ggplot(chemo_multi_melt, aes(x=Buoyant_density, y=rel, fill=Isotope, color=Isotope)) +
  geom_point() +
  geom_line() +
  theme_classic2(base_size = 10) +
  facet_wrap(~Strain) +
 # scale_y_log10() + 
  scale_fill_manual(values=colors_fill) +
  scale_color_manual(values=colors_fill) +
  ylab("Ratio of Maximum Quantity") +
  xlab("Buoyant Density (g/mL)") +
  theme(legend.position = 'bottom')


