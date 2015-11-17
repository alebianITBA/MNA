function frecuency_from_mel = mel_inverse(value)
	frecuency_from_mel = 700* (e**(value/1125)-1);
end