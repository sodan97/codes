const dateFormatter = new Intl.DateTimeFormat("fr-FR", { day: "numeric", month: "long", year: "numeric", timeZone: "UTC" });

export const formatDate = (iso: string) => dateFormatter.format(new Date(`${iso}T00:00:00Z`));
