#!/usr/bin/python3
import math

# Calculate time on air using various SFs
# Suppose the payload is 50 Bytes
for sf in [7, 8, 9, 10]:
	nPreamble = 8
	bwHz = 125000 # Hz
	pl = 50 # Bytes
	cr = 1 # Coding rate (1 corresponding to 4/5, 4 to 4/8)
	ih = 0	# Whether header is enabled
	crc = 1 # Whether CRC is enabled
	de = 0 # Low data rate optimize

	tSym = pow(2, sf) / bwHz
	tPreamble = (nPreamble + 4.25) * tSym

	num = 8 * pl - 4 * sf + 28 + 16 * crc - 20 * ih;
	den = 4 * (sf - 2 * de)
	payloadSymbNb = 8 + max (math.ceil (num / den) * (cr + 4), .0);
	tPayload = payloadSymbNb * tSym

	print('SF{}, # of payload symbols: {}'.format(sf, payloadSymbNb))
	print('Preamble time: {}, payload time: {}, total time-on-air: {}'.format(\
		tPreamble, tPayload, tPreamble + tPayload))
