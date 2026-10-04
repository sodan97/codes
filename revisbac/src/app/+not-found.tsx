import { Redirect } from 'expo-router';

/** Adresse inconnue (lien périmé, ancienne version) : retour à l'accueil. */
export default function NotFound() {
  return <Redirect href="/" />;
}
