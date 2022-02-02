import telebot

#add your telegram bot API key here
API_KEY = ""


bot = telebot.TeleBot(API_KEY)
DEREFERER_BASE_URL = "https://dereferer.me/?"


def isURL(message):
    return message.text.startswith('https') and not message.text.startswith(DEREFERER_BASE_URL)
    

@bot.message_handler(func=isURL)
def dereferer(message):
    url = DEREFERER_BASE_URL + message.text.replace(':', '%3A')
    bot.send_message(message.chat.id, url)
    bot.delete_message(message.chat.id, message.message_id)

bot.infinity_polling()
