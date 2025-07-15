raw_matrix = dat.X.todense()
scrub = scr.Scrublet(raw_matrix)
doublet_score, predicted_doublets = scrub.scrub_doublets()
real_cells = np.logical_not(predicted_doublets)
dat = dat[real_cells,:]
